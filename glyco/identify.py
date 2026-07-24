"""
동정 + 결과 조립
----------------
흐름:
  1) 후보 조성 생성 -> adduct 타깃 m/z 목록
  2) MS2 precursor 를 타깃에 ppm 매칭 -> '동정된' 글리칸 확정(=MS/MS 근거 있음)
  3) 동정된 글리칸의 각 adduct XIC apex 합산 -> 정량
  4) 상대량(%) + 구조 분류

설계 결정(advisor):
  - 동정은 결정론적 m/z 매칭으로 확신 있게.
  - 정량(%)은 사람이 읽은 NL 값과 절대값이 다르므로 순위/패턴으로 해석.
  - RT 는 시트값을 쓰지 않고 이번 run 의 XIC 에서 검출.
"""

import numpy as np

from . import masses, compositions, quantify, classify, diagnostic
from .chem import ppm_error


def match_ms2(msdata, mz_targets, meta, ppm_tol=10.0):
    """
    각 MS2 precursor 를 가장 가까운 타깃에 매칭.
    반환: dict[cand_idx] -> list of dict(scan, rt, adduct, z, obs_mz, ppm)
    """
    from collections import defaultdict
    hits = defaultdict(list)
    if len(mz_targets) == 0:
        return hits
    for rt, pmz, z, scan in msdata.ms2:
        idx = np.searchsorted(mz_targets, pmz)
        cands = []
        if idx < len(mz_targets):
            cands.append(idx)
        if idx > 0:
            cands.append(idx - 1)
        best = None
        for j in cands:
            ppm = abs(ppm_error(pmz, mz_targets[j]))
            if ppm <= ppm_tol and (best is None or ppm < best[0]):
                best = (ppm, j)
        if best is None:
            continue
        j = best[1]
        m = meta[j]
        # MS2 의 보고된 charge 가 있으면 일치하는 것만 신뢰
        if z is not None and z != m["z"]:
            continue
        hits[m["cand_idx"]].append({
            "scan": scan, "rt": rt, "adduct": m["adduct"], "z": m["z"],
            "obs_mz": pmz, "ppm": ppm_error(pmz, m["mz"]),
        })
    return hits


def screening_table(msdata, diag_table, ppm=20.0, anchor="HexNAc"):
    """
    Xcalibur 수작업 '스크리닝 시트' 재현.
    각 MS2 스캔에서 진단이온(oxonium)들의 관측 m/z + precursor/charge/RT 를 뽑는다.
    anchor 이온(기본 HexNAc 204)이 보이는 스캔(=글리칸 후보)만 반환.
    anchor 는 문자열 또는 문자열 리스트이며, 리스트면 하나라도 보이는 스캔을 반환한다.
    반환: (rows, ion_names)  — rows[i] = dict(scan, rt, precursor, charge, ions={name:obs_mz})
    """
    import numpy as np
    ion_names = [n for n, _ in diag_table]
    if isinstance(anchor, str):
        anchors = [x.strip() for x in anchor.split(",") if x.strip()]
    else:
        anchors = list(anchor or [])
    rows = []
    for rt, pmz, z, scan in msdata.ms2:
        pk = msdata.ms2_peaks.get(scan)
        if pk is None:
            continue
        mz = np.sort(pk[0])
        ions = {name: diagnostic.nearest_within_ppm(mz, target, ppm) for name, target in diag_table}
        if anchors and not any(ions.get(a) is not None for a in anchors):
            continue
        info = getattr(msdata, "ms2_info", {}).get(scan, {})
        rows.append({
            "scan": scan,
            "rt": rt,
            "precursor": pmz,
            "isolation_target": info.get("isolation_target_mz"),
            "monoisotope": info.get("selected_mz", pmz),
            "charge": z,
            "ions": ions,
        })
    rows.sort(key=lambda r: r["rt"])
    return rows, ion_names


def screen_diagnostics(msdata, spec):
    """
    진단규칙(spec)으로 MS2 스캔을 스크리닝(단일 글리칸-공통 규칙 + feature 주석).
    - 코어이온(spec.ions) 관측 m/z·ppm오차 추출.
    - 채택 = spec.accept 의 모든 그룹 만족(그룹=any 이온 중 min개 이상 검출).
    - feature = spec.features 각각(리스트면 OR)이 검출되면 이름 주석.
    반환: rows = [{scan, rt, precursor, monoisotope, charge,
                   ion_obs={name:(obs|None, ppm|None)}, features=[name,...]}]  (채택분만)
    """
    import numpy as np
    rows = []
    for rt, pmz, z, scan in msdata.ms2:
        pk = msdata.ms2_peaks.get(scan)
        if pk is None:
            continue
        mz = np.sort(pk[0])
        obs = {}
        for name, theo in spec.ions.items():
            o = diagnostic.nearest_within_ppm(mz, theo, spec.ppm)
            obs[name] = (o, ppm_error(o, theo))
        # 채택 판정: 모든 그룹에서 (검출된 any 이온 수) >= min
        accepted = all(
            sum(1 for n in grp["any"] if obs.get(n, (None,))[0] is not None) >= grp["min"]
            for grp in spec.accept
        )
        if not accepted:
            continue
        feats = [ft["name"] for ft in spec.features
                 if any(diagnostic.nearest_within_ppm(mz, m, spec.ppm) is not None for m in ft["mz"])]
        info = getattr(msdata, "ms2_info", {}).get(scan, {})
        rows.append({
            "scan": scan, "rt": rt, "precursor": pmz,
            "monoisotope": info.get("monoisotopic_mz") or pmz,
            "charge": z, "ion_obs": obs, "features": feats,
        })
    rows.sort(key=lambda r: r["rt"])
    for i, r in enumerate(rows):
        r["sc_no"] = i + 1          # Screening 행 번호(구조찾기에서 역참조)
    return rows


def run(candidates, *, msdata, max_charge=3, ppm_tol=10.0,
        ms1_ppm=5.0, require_ms2=True, min_intensity=0.0,
        quant_method="area", rt_window=0.5, rt_consistency=1.0,
        chem=None, require_diagnostic=None, ms2_ppm=20.0,
        ms1_first=False, ms1_min_adducts=3, ms1_first_ppm=3.0,
        ms1_noise_factor=70.0, log=print):
    """전체 동정+정량 파이프라인. 반환: 결과 dict 리스트(상대량 내림차순).

    quant_method: 'area'(EIC 면적, 논문방식) | 'apex'(피크높이)
    chem + require_diagnostic: 진단 oxonium 확인 활성화(MS2 피크 보관 필요).
    ms1_first: MS2 없이도 '같은 RT에 co-eluting adduct ≥ ms1_min_adducts' 이면
               글리칸으로 인정(시알산 등 조각화 안 된 것 회수). MS2 경로와 합집합.
    """
    adducts = compositions.adduct_ions(max_charge=max_charge)
    ion_mz_fn = chem.ion_mz if chem is not None else masses.ion_mz
    mz_targets, meta = quantify.build_targets(candidates, adducts, ion_mz_fn)
    log(f"[동정] 후보 {len(candidates)}개 × adduct -> 타깃 {len(mz_targets)}개")

    ms2_hits = match_ms2(msdata, mz_targets, meta, ppm_tol=ppm_tol)
    log(f"[동정] MS2 근거가 있는 후보 {len(ms2_hits)}개")

    # --- 경로 1 (MS2): precursor 매칭 + (선택) 진단 oxonium 확인 ---
    if require_ms2:
        ms2_keep = set(ms2_hits.keys())
    else:
        ms2_keep = set(range(len(candidates)))
    if require_diagnostic and chem is not None and msdata.ms2_peaks:
        diag_table = chem.diagnostic_table()
        before = len(ms2_keep)
        confirmed = set()
        for ci in list(ms2_keep):
            for h in ms2_hits.get(ci, []):
                if diagnostic.confirm_scan(msdata, h["scan"], diag_table,
                                           require_diagnostic, ppm_tol=ms2_ppm):
                    confirmed.add(ci)
                    break
        ms2_keep = confirmed
        log(f"[동정] 진단이온 확인 통과 {len(ms2_keep)}/{before}개")

    # --- 경로 2 (MS1-first): co-eluting adduct 다수로 인정 (MS2 불필요) ---
    ms1_keep = set()
    ms1_evidence = {}
    if ms1_first:
        # 강도 하한 = 노이즈 median × factor.
        # (노이즈는 장비 특성이라 시료 로딩·동적범위와 무관 → 시료 간 일반화됨)
        import numpy as _np
        allint = _np.concatenate([it for _, _, it in msdata.ms1 if it.size]) \
            if msdata.ms1 else _np.array([0.0])
        noise = float(_np.median(allint)) if allint.size else 0.0
        floor = ms1_noise_factor * noise
        a_int, a_rt, a_mz = quantify.apex_scan(msdata, mz_targets, ppm_tol=ms1_first_ppm, log=log)
        from collections import defaultdict
        by_cand = defaultdict(list)   # ci -> [(intensity, rt), ...] (정밀질량·강도 통과분)
        for i, m in enumerate(meta):
            amz = a_mz[i]
            if amz != amz or a_int[i] <= floor:
                continue
            if abs(ppm_error(amz, m["mz"])) > ms1_first_ppm:
                continue
            by_cand[m["cand_idx"]].append((a_int[i], a_rt[i]))
        for ci, lst in by_cand.items():
            dom_rt = max(lst, key=lambda x: x[0])[1]     # 최강 adduct 의 RT
            n_co = sum(1 for _, rt in lst if abs(rt - dom_rt) <= rt_consistency)
            if n_co >= ms1_min_adducts:
                ms1_keep.add(ci)
                ms1_evidence[ci] = n_co
        log(f"[동정] MS1-first (co-eluting adduct≥{ms1_min_adducts}) 통과 {len(ms1_keep)}개")

    keep = ms2_keep | ms1_keep
    if not keep:
        log("[경고] 동정된 글리칸이 없습니다 (ppm 허용오차/범위를 확인하세요).")
        return []

    # 타깃 중 keep 에 속한 것만 골라 XIC (동정 후 소수 타깃만 면적적분)
    sel = [i for i, m in enumerate(meta) if m["cand_idx"] in keep]
    sel_mz = mz_targets[sel]
    xic = quantify.extract_xic(msdata, sel_mz, ppm_tol=ppm_tol, log=log)
    if quant_method == "apex":
        signal, apex_rt, apex_mz = quantify.apex(xic)
    else:
        signal, apex_rt, apex_mz = quantify.area(xic, rt_window=rt_window)

    # 1차: cand_idx -> 검출된 adduct 후보 수집(아직 합산 전)
    from collections import defaultdict
    raw = defaultdict(list)
    for local_i, target_i in enumerate(sel):
        m = meta[target_i]
        inten = float(signal[local_i])
        if inten <= min_intensity:
            continue
        # MS1 정밀질량 게이트: apex 피크가 이론 m/z 와 ms1_ppm 이내여야 인정
        amz = float(apex_mz[local_i])
        if amz != amz or abs(ppm_error(amz, m["mz"])) > ms1_ppm:   # NaN 또는 ppm 초과
            continue
        raw[m["cand_idx"]].append({
            "adduct": compositions.adduct_label(m["adduct"]), "z": m["z"],
            "mz": m["mz"], "intensity": inten, "rt": float(apex_rt[local_i]),
        })

    # 2차: RT 일치 필터 — 주 피크(최대강도)와 같은 RT 의 adduct만 합산
    #      (같은 m/z 의 다른 co-eluting 분자를 잘못 주워담는 것 방지)
    per = {}
    for ci, ads in raw.items():
        main = max(ads, key=lambda a: a["intensity"])
        if rt_consistency:
            kept = [a for a in ads if abs(a["rt"] - main["rt"]) <= rt_consistency]
        else:
            kept = ads
        per[ci] = {"sum": sum(a["intensity"] for a in kept), "adducts": kept,
                   "rt": main["rt"], "best": main["intensity"]}

    total = sum(v["sum"] for v in per.values())
    results = []
    for ci, v in per.items():
        if v["sum"] <= 0:
            continue
        cand = candidates[ci]
        c = cand["composition"]
        cls = classify.classify(c)
        m2 = ms2_hits.get(ci, [])
        ev = []
        if ci in ms2_keep:
            ev.append("MS2")
        if ci in ms1_keep:
            ev.append("MS1")
        best_adduct = max(v["adducts"], key=lambda a: a["intensity"]) if v["adducts"] else None
        results.append({
            "name": cand["name"],
            "oxford": cand.get("oxford", ""),
            "composition": c,
            "formula": cand["formula"],
            "neutral": cand["neutral"],
            "type": cls["type"],
            "sialylated": cls["sialylated"],
            "fucosylated": cls["fucosylated"],
            "n_sialic": cls["n_sialic"],
            "n_fucose": cls["n_fucose"],
            "best_mz": best_adduct["mz"] if best_adduct else None,
            "best_adduct": best_adduct["adduct"] if best_adduct else None,
            "best_z": best_adduct["z"] if best_adduct else None,
            "rt": v["rt"],
            "ms2_count": len(m2),
            "evidence": "+".join(ev) or "MS1",
            "ms1_adducts": ms1_evidence.get(ci),
            "best_ppm": min((abs(h["ppm"]) for h in m2), default=None),
            "intensity_sum": v["sum"],
            "relative_pct": v["sum"] / total * 100 if total else 0.0,
            "adducts": v["adducts"],
        })
    results.sort(key=lambda r: r["relative_pct"], reverse=True)
    log(f"[동정] 최종 정량 글리칸 {len(results)}개")
    return results
