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
    p.add_argument("--no-ms2", action="store_true", help="MS2 근거 없이도 정량")
    p.add_argument("--no-diagnostic", action="store_true", help="진단 oxonium 확인 끄기")
    p.add_argument("--min-intensity", type=float, help="이 값 미만 adduct 무시")
    p.add_argument("--keep-mzml", action="store_true", help="변환 mzML 보존")
    args = p.parse_args(argv)

    if not os.path.isfile(args.input):
        print(f"[오류] 입력 파일이 없습니다: {args.input}")
        return 2

    overrides = {}
    if args.precursor_ppm is not None: overrides["precursor_ppm"] = args.precursor_ppm
    if args.ms1_ppm is not None: overrides["ms1_ppm"] = args.ms1_ppm
    if args.max_charge is not None: overrides["max_charge"] = args.max_charge
    if args.quant: overrides["quant_method"] = args.quant
    if args.rt_window is not None: overrides["rt_window"] = args.rt_window
    if args.no_ms2: overrides["require_ms2"] = False
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
