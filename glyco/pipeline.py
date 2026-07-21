"""
파이프라인 오케스트레이션
-------------------------
설정 + 입력 파일 -> 변환 -> 파싱 -> 조성생성 -> 동정 -> 정량 -> 결과.
CLI/GUI/노트북 어디서든 이 run() 하나만 호출하면 된다.
"""

import glob
import os

from . import config as config_mod
from . import raw, mzml_parse, compositions, identify, report, aggregate, targets as targets_mod
from .chem import Chemistry


def targeted_screening(input_path, targets_file, *, config_path=None, output=None,
                       work_dir=None, keep_mzml=False, overrides=None, log=print):
    """
    사용자 지정 진단이온 목록(targets_file)으로 MS2 타깃 스크리닝.
    반환: (rows, out_path, spec)
    """
    overrides = overrides or {}
    spec = targets_mod.load(targets_file)
    if "ppm" in overrides: spec.ppm = overrides["ppm"]
    if "min_hits" in overrides: spec.min_hits = overrides["min_hits"]
    if "precursor_floor" in overrides: spec.precursor_floor = overrides["precursor_floor"]
    targets_mod._validate(spec)
    log(targets_mod.summary(spec))

    out = output or (os.path.splitext(input_path)[0] + "_targeted.xlsx")
    work = work_dir or os.path.join(os.path.dirname(os.path.abspath(input_path)), "_glycan_targeted_work")
    mzml_path = raw.to_mzml(input_path, work, ms_levels="2", log=log)
    data = mzml_parse.parse(mzml_path, keep_ms2_peaks=True, keep_ms1=False, log=log)

    rows = identify.targeted_screen(data, spec.glycans, ppm=spec.ppm, min_hits=spec.min_hits)
    nm = sum(1 for r in rows if r["tier"] == "matched")
    nh = sum(1 for r in rows if r["tier"] == "holding")
    log(f"[타깃] 매칭(≥{spec.min_hits}) {nm} | 보류(부분) {nh} 쌍")

    report.write_targeted(rows, spec, out, meta={"sample": os.path.basename(input_path)})
    if not keep_mzml and mzml_path.startswith(work):
        try:
            os.remove(mzml_path); os.rmdir(work)
        except OSError:
            pass
    return rows, out, spec


def find_inputs(directory):
    """디렉토리에서 .raw / .mzML 입력을 찾아 정렬 반환(재귀 X)."""
    fs = []
    for ext in ("*.raw", "*.mzML", "*.mzml"):
        fs += glob.glob(os.path.join(directory, ext))
    return sorted(set(fs))


def batch(input_dir, *, config_path=None, output_dir=None, overrides=None,
          keep_mzml=False, log=print):
    """
    디렉토리의 여러 .raw/.mzML 을 각각 분석(개별 엑셀)하고, 반복으로 취합한다.
    반환: (per_file[(label,results)], agg, agg_path)
    """
    files = find_inputs(input_dir)
    if not files:
        raise FileNotFoundError(f"{input_dir} 에 .raw/.mzML 이 없습니다.")
    out_dir = output_dir or os.path.join(os.path.abspath(input_dir), "glycoquant_results")
    os.makedirs(out_dir, exist_ok=True)

    per_file = []
    failed = []
    for i, f in enumerate(files, 1):
        label = os.path.splitext(os.path.basename(f))[0]
        log(f"\n===== [{i}/{len(files)}] {label} =====")
        per_out = os.path.join(out_dir, label + "_glycans.xlsx")
        try:
            results, _, cfg = analyze(f, config_path=config_path, output=per_out,
                                      keep_mzml=keep_mzml, overrides=overrides, log=log)
            per_file.append((label, results or []))
        except Exception as e:   # 파일 하나 실패해도 나머지는 계속
            failed.append((label, str(e)))
            log(f"[경고] '{label}' 처리 실패 — 건너뜀: {e}")

    if not per_file:
        raise RuntimeError("모든 파일 처리 실패: " + "; ".join(f"{n}({e})" for n, e in failed))
    if failed:
        log(f"[주의] {len(failed)}개 파일 실패(취합서 제외): " + ", ".join(n for n, _ in failed))

    log(f"\n===== 반복 취합 (n={len(per_file)}) =====")
    agg = aggregate.aggregate(per_file)
    agg_path = os.path.join(out_dir, "aggregated.xlsx")
    report.write_aggregated(agg, agg_path, meta={"group": os.path.basename(os.path.abspath(input_dir))})
    log(f"[취합] 글리칸 {len(agg['glycans'])}종 -> {agg_path}")
    return per_file, agg, agg_path


def analyze(input_path, *, config_path=None, output=None, work_dir=None,
            keep_mzml=False, overrides=None, log=print):
    """
    input_path : .raw 또는 .mzML
    config_path: YAML 설정(없으면 기본 ProA)
    overrides  : 설정을 덮어쓸 dict (예: {'ms1_ppm':3, 'ranges':{...}})
    반환: (results, output_path, cfg)
    """
    overrides = overrides or {}
    cfg = config_mod.load(config_path)
    chem = Chemistry(cfg)

    out = output or (os.path.splitext(input_path)[0] + "_glycans.xlsx")
    work = work_dir or os.path.join(os.path.dirname(os.path.abspath(input_path)), "_glycan_work")

    # 1) raw -> mzML
    mzml_path = raw.to_mzml(input_path, work, log=log)

    # 2) 파싱 (진단이온 확인 또는 스크리닝 시트 → MS2 피크 보관)
    use_diag = bool(cfg.require_diagnostic) and not overrides.get("no_diagnostic")
    want_screening = not overrides.get("no_screening")
    data = mzml_parse.parse(mzml_path, keep_ms2_peaks=(use_diag or want_screening), log=log)

    # 3) 조성 생성 (설정의 라벨 화학으로 질량 계산)
    if "ranges" in overrides and overrides["ranges"]:
        cands = compositions.generate(ranges=overrides["ranges"], rules=cfg.plausibility,
                                      label_count=cfg.label.count, chem=chem)
    else:
        cands = compositions.from_config(cfg, chem=chem)
    log(f"[조성] 후보 {len(cands)}개 생성")

    # 4) 동정 + 정량
    q = cfg.quantify
    mf = cfg.ms1_first
    ms1_first_on = overrides.get("ms1_first", mf.get("enabled", False))
    results = identify.run(
        cands, msdata=data,
        max_charge=overrides.get("max_charge", cfg.max_charge),
        ppm_tol=overrides.get("precursor_ppm", cfg.precursor_ppm),
        ms1_ppm=overrides.get("ms1_ppm", cfg.ms1_ppm),
        require_ms2=overrides.get("require_ms2", q.get("require_ms2", True)),
        min_intensity=overrides.get("min_intensity", 0.0),
        quant_method=overrides.get("quant_method", q.get("method", "area")),
        rt_window=overrides.get("rt_window", q.get("rt_window_min", 0.5)),
        rt_consistency=overrides.get("rt_consistency", q.get("rt_consistency_min", 1.0)),
        chem=chem,                                  # 질량+진단 모두 설정 화학 사용
        require_diagnostic=cfg.require_diagnostic if use_diag else None,
        ms2_ppm=cfg.ms2_ppm,
        ms1_first=ms1_first_on,
        ms1_min_adducts=overrides.get("ms1_min_adducts", mf.get("min_adducts", 3)),
        ms1_first_ppm=overrides.get("ms1_first_ppm", mf.get("ppm", 3.0)),
        ms1_noise_factor=overrides.get("ms1_noise_factor", mf.get("noise_factor", 70.0)),
        log=log,
    )

    # 5) 스크리닝 시트(Xcalibur 수작업 대체) + 리포트
    screening = screening_ions = None
    if want_screening and data.ms2_peaks:
        # screening_all=True 면 필터 없이 모든 MS2 스캔을 그대로 덤프
        anchor = None if overrides.get("screening_all") else "HexNAc"
        screening, screening_ions = identify.screening_table(
            data, chem.diagnostic_table(), anchor=anchor)
        kind = "전체 MS2" if anchor is None else "글리칸(204+)"
        log(f"[스크리닝] {kind} 스캔 {len(screening)}개 표로 추출")
    if results:
        meta = {"sample": os.path.basename(input_path), "source": cfg.name}
        report.write(results, out, meta=meta, screening=screening, screening_ions=screening_ions)

    # 정리
    if not keep_mzml and mzml_path.startswith(work):
        try:
            os.remove(mzml_path); os.rmdir(work)
        except OSError:
            pass

    return results, out, cfg


def screening(input_path, *, config_path=None, output=None, work_dir=None,
              keep_mzml=False, anchor="HexNAc,ProA-HexNAc", ppm=None, log=print):
    """
    .raw 또는 .mzML 에서 MS2 진단이온 스크리닝 시트만 생성한다.

    Xcalibur 에서 scan/RT/oxonium/precursor 값을 손으로 옮기던 작업을
    대체하기 위한 빠른 경로다. .raw 입력은 MS2만 mzML 로 변환한다.
    """
    cfg = config_mod.load(config_path)
    chem = Chemistry(cfg)

    out = output or (os.path.splitext(input_path)[0] + "_screening.xlsx")
    work = work_dir or os.path.join(os.path.dirname(os.path.abspath(input_path)), "_glycan_screening_work")

    mzml_path = raw.to_mzml(input_path, work, ms_levels="2", log=log)
    data = mzml_parse.parse(mzml_path, keep_ms2_peaks=True, keep_ms1=False, log=log)

    screening_rows, screening_ions = identify.screening_table(
        data, chem.diagnostic_table(), ppm=ppm or cfg.ms2_ppm, anchor=anchor)
    log(f"[스크리닝] 글리칸 후보 MS2 스캔 {len(screening_rows)}개 추출")

    meta = {"sample": os.path.basename(input_path), "source": cfg.name}
    report.write_screening(screening_rows, screening_ions, out, meta=meta)

    if not keep_mzml and mzml_path.startswith(work):
        try:
            os.remove(mzml_path); os.rmdir(work)
        except OSError:
            pass

    return screening_rows, out, cfg
