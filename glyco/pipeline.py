"""
파이프라인 오케스트레이션
-------------------------
설정 + 입력 파일 -> 변환 -> 파싱 -> 조성생성 -> 동정 -> 정량 -> 결과.
CLI/GUI/노트북 어디서든 이 run() 하나만 호출하면 된다.
"""

import os

from . import config as config_mod
from . import raw, mzml_parse, compositions, identify, report
from .chem import Chemistry


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

    # 2) 파싱 (진단이온 확인 시 MS2 피크 보관)
    use_diag = bool(cfg.require_diagnostic) and not overrides.get("no_diagnostic")
    data = mzml_parse.parse(mzml_path, keep_ms2_peaks=use_diag, log=log)

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

    # 5) 리포트
    if results:
        meta = {"sample": os.path.basename(input_path), "source": cfg.name}
        report.write(results, out, meta=meta)

    # 정리
    if not keep_mzml and mzml_path.startswith(work):
        try:
            os.remove(mzml_path); os.rmdir(work)
        except OSError:
            pass

    return results, out, cfg
