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

from . import masses, compositions, quantify, classify


def _ppm(obs, theo):
    return (obs - theo) / theo * 1e6


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
            ppm = abs(_ppm(pmz, mz_targets[j]))
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
            "obs_mz": pmz, "ppm": _ppm(pmz, m["mz"]),
        })
    return hits


def run(candidates, *, msdata, max_charge=3, ppm_tol=10.0,
        ms1_ppm=5.0, require_ms2=True, min_intensity=0.0, log=print):
    """전체 동정+정량 파이프라인. 반환: 결과 dict 리스트(상대량 내림차순)."""
    adducts = compositions.adduct_ions(max_charge=max_charge)
    mz_targets, meta = quantify.build_targets(candidates, adducts, masses.ion_mz)
    log(f"[동정] 후보 {len(candidates)}개 × adduct -> 타깃 {len(mz_targets)}개")

    ms2_hits = match_ms2(msdata, mz_targets, meta, ppm_tol=ppm_tol)
    log(f"[동정] MS2 근거가 있는 후보 {len(ms2_hits)}개")

    # 정량 대상 결정
    if require_ms2:
        keep = set(ms2_hits.keys())
    else:
        keep = set(range(len(candidates)))
    if not keep:
        log("[경고] 동정된 글리칸이 없습니다 (ppm 허용오차/범위를 확인하세요).")
        return []

    # 타깃 중 keep 에 속한 것만 골라 XIC
    sel = [i for i, m in enumerate(meta) if m["cand_idx"] in keep]
    sel_mz = mz_targets[sel]
    apex_int, apex_rt, apex_mz = quantify.xic_apex(msdata, sel_mz, ppm_tol=ppm_tol, log=log)

    # cand_idx -> adduct별 강도 집계
    from collections import defaultdict
    per = defaultdict(lambda: {"sum": 0.0, "adducts": [], "rt": None, "best": 0.0})
    for local_i, target_i in enumerate(sel):
        m = meta[target_i]
        inten = float(apex_int[local_i])
        if inten <= min_intensity:
            continue
        # MS1 정밀질량 게이트: apex 피크가 이론 m/z 와 ms1_ppm 이내여야 인정
        amz = float(apex_mz[local_i])
        if amz != amz or abs(_ppm(amz, m["mz"])) > ms1_ppm:   # NaN 또는 ppm 초과
            continue
        ci = m["cand_idx"]
        per[ci]["sum"] += inten
        per[ci]["adducts"].append({
            "adduct": compositions.adduct_label(m["adduct"]), "z": m["z"],
            "mz": m["mz"], "intensity": inten, "rt": float(apex_rt[local_i]),
        })
        if inten > per[ci]["best"]:
            per[ci]["best"] = inten
            per[ci]["rt"] = float(apex_rt[local_i])

    total = sum(v["sum"] for v in per.values())
    results = []
    for ci, v in per.items():
        if v["sum"] <= 0:
            continue
        cand = candidates[ci]
        c = cand["composition"]
        cls = classify.classify(c)
        m2 = ms2_hits.get(ci, [])
        best_adduct = max(v["adducts"], key=lambda a: a["intensity"]) if v["adducts"] else None
        results.append({
            "name": cand["name"],
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
            "best_ppm": min((abs(h["ppm"]) for h in m2), default=None),
            "intensity_sum": v["sum"],
            "relative_pct": v["sum"] / total * 100 if total else 0.0,
            "adducts": v["adducts"],
        })
    results.sort(key=lambda r: r["relative_pct"], reverse=True)
    log(f"[동정] 최종 정량 글리칸 {len(results)}개")
    return results
