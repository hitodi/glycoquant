#!/usr/bin/env python3
"""
glycan_analyze — 라벨 글리칸 분석기 (Thermo .raw / .mzML -> 동정+정량 Excel)

기본은 ProA N-글리칸 설정(configs/proa_nglycan.yaml).
다른 라벨/단당/범위는 --config 로 YAML 만 바꿔 끼우면 됨.

사용법:
    python glycan_analyze.py 시료.raw
    python glycan_analyze.py 시료.raw -o 결과.xlsx --ms1-ppm 3
    python glycan_analyze.py 시료.raw --config configs/2ab_nglycan.yaml
"""
import argparse
import os
import sys
import time

from glyco import pipeline


def main(argv=None):
    p = argparse.ArgumentParser(
        description="라벨 글리칸 분석 (Thermo .raw/.mzML -> 동정+정량 Excel)")
    p.add_argument("input", help="입력 .raw 또는 .mzML")
    p.add_argument("-o", "--output", help="출력 .xlsx (기본: 입력명_glycans.xlsx)")
    p.add_argument("-c", "--config", help="설정 YAML (기본: configs/proa_nglycan.yaml)")
    p.add_argument("--precursor-ppm", type=float, help="MS2 precursor 매칭 ppm")
    p.add_argument("--ms1-ppm", type=float, help="MS1 정밀질량 게이트 ppm")
    p.add_argument("--max-charge", type=int, help="최대 전하수")
    p.add_argument("--quant", choices=["area", "apex"], help="정량 방식")
    p.add_argument("--rt-window", type=float, help="면적적분 RT 윈도우(분)")
    p.add_argument("--rt-consistency", type=float, help="adduct 합산 RT 일치 허용(분); 0=끔")
    p.add_argument("--ms1-first", action="store_true",
                   help="MS1-first 모드: 조각화 안 된(시알산 등) 글리칸을 co-eluting adduct로 회수. "
                        "⚠️정밀도↓/회수↑ 트레이드오프 — 시알산 풍부 시료에서만 권장")
    p.add_argument("--ms1-min-adducts", type=int, help="MS1-first: 같은 RT 공존 adduct 최소 수(기본 3)")
    p.add_argument("--ms1-first-ppm", type=float, help="MS1-first: MS1 정밀질량 ppm(기본 3)")
    p.add_argument("--ms1-noise-factor", type=float, help="MS1-first: 강도하한=노이즈median×이값(기본 70)")
    p.add_argument("--no-ms2", action="store_true", help="MS2 근거 없이도 정량")
    p.add_argument("--no-diagnostic", action="store_true", help="진단 oxonium 확인 끄기")
    p.add_argument("--min-intensity", type=float, help="이 값 미만 adduct 무시")
    p.add_argument("--no-screening", action="store_true",
                   help="스크리닝 시트(스캔별 진단이온/precursor) 생략")
    p.add_argument("--screening-all", action="store_true",
                   help="스크리닝을 필터 없이 전체 MS2 스캔으로(204 없는 스캔도 전부 덤프)")
    p.add_argument("--screening-only", "--screen-only", action="store_true",
                   help="동정·정량 없이 .raw/.mzML에서 Screening 엑셀만 생성")
    p.add_argument("--screening-anchor", default="HexNAc,ProA-HexNAc",
                   help="Screening에 포함할 anchor 진단이온 이름들(쉼표 구분, 기본: HexNAc,ProA-HexNAc; 빈 문자열이면 전체 MS2)")
    p.add_argument("--screening-ppm", type=float,
                   help="Screening 진단이온 m/z 허용오차 ppm(기본: 설정 ms2_ppm)")
    p.add_argument("--keep-mzml", action="store_true", help="변환 mzML 보존")
    args = p.parse_args(argv)

    if not os.path.isfile(args.input):
        print(f"[오류] 입력 파일이 없습니다: {args.input}")
        return 2

    if args.screening_only:
        t0 = time.time()
        anchor = None if args.screening_all else (args.screening_anchor or None)
        rows, out, cfg = pipeline.screening(
            args.input, config_path=args.config, output=args.output,
            keep_mzml=args.keep_mzml, anchor=anchor, ppm=args.screening_ppm)
        print(f"\n[완료] 설정='{cfg.name}'  Screening 스캔 {len(rows)}개 -> {out}  ({time.time()-t0:.1f}s)")
        return 0

    overrides = {}
    if args.precursor_ppm is not None: overrides["precursor_ppm"] = args.precursor_ppm
    if args.ms1_ppm is not None: overrides["ms1_ppm"] = args.ms1_ppm
    if args.max_charge is not None: overrides["max_charge"] = args.max_charge
    if args.quant: overrides["quant_method"] = args.quant
    if args.rt_window is not None: overrides["rt_window"] = args.rt_window
    if args.rt_consistency is not None:
        overrides["rt_consistency"] = args.rt_consistency or None
    if args.ms1_first: overrides["ms1_first"] = True
    if args.ms1_min_adducts is not None: overrides["ms1_min_adducts"] = args.ms1_min_adducts
    if args.ms1_first_ppm is not None: overrides["ms1_first_ppm"] = args.ms1_first_ppm
    if args.ms1_noise_factor is not None: overrides["ms1_noise_factor"] = args.ms1_noise_factor
    if args.no_ms2: overrides["require_ms2"] = False
    if args.no_screening: overrides["no_screening"] = True
    if args.screening_all: overrides["screening_all"] = True
    if args.no_diagnostic: overrides["no_diagnostic"] = True
    if args.min_intensity is not None: overrides["min_intensity"] = args.min_intensity

    t0 = time.time()
    results, out, cfg = pipeline.analyze(
        args.input, config_path=args.config, output=args.output,
        keep_mzml=args.keep_mzml, overrides=overrides)

    if not results:
        print("[결과] 동정된 글리칸이 없습니다. --ms1-ppm 를 늘리거나 --no-diagnostic/--no-ms2 를 시도하세요.")
        return 1

    print(f"\n[완료] 설정='{cfg.name}'  글리칸 {len(results)}개 -> {out}  ({time.time()-t0:.1f}s)")
    print("상위 5개:")
    for r in results[:5]:
        print(f"  {r['oxford']:12s} {r['name']:30s} {r['type']:12s} "
              f"{r['relative_pct']:5.2f}%  (m/z {r['best_mz']:.4f} {r['best_adduct']} {r['best_z']}+)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
