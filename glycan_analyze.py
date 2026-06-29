#!/usr/bin/env python3
"""
glycan_analyze — ProA 표지 N-글리칸 분석기 (Thermo .raw -> Excel)

사용법:
    python glycan_analyze.py 시료.raw
    python glycan_analyze.py 시료.raw -o 결과.xlsx --ppm 8 --max-charge 3

비개발자용: run.command(맥) / run.bat(윈도우) 더블클릭 후 .raw 경로 입력.
"""

import argparse
import os
import sys
import time

from glyco import raw, mzml_parse, compositions, identify, report


def build_ranges(args):
    r = compositions.SearchRanges()
    if args.hexnac:
        r.HexNAc = range(args.hexnac[0], args.hexnac[1] + 1)
    if args.hex:
        r.Hex = range(args.hex[0], args.hex[1] + 1)
    if args.fuc:
        r.dHex = range(args.fuc[0], args.fuc[1] + 1)
    if args.neu5ac:
        r.Neu5Ac = range(args.neu5ac[0], args.neu5ac[1] + 1)
    if args.neu5gc:
        r.Neu5Gc = range(args.neu5gc[0], args.neu5gc[1] + 1)
    if args.xyl:
        r.Xyl = range(0, 2)
    return r


def main(argv=None):
    p = argparse.ArgumentParser(
        description="ProA 표지 N-글리칸 분석 (Thermo .raw / .mzML -> 동정+정량 Excel)")
    p.add_argument("input", help="입력 .raw 또는 .mzML 파일")
    p.add_argument("-o", "--output", help="출력 .xlsx (기본: 입력명_glycans.xlsx)")
    p.add_argument("--ppm", type=float, default=10.0,
                   help="MS2 precursor 매칭 허용오차 ppm (isolation 값이라 느슨, 기본 10)")
    p.add_argument("--ms1-ppm", type=float, default=5.0,
                   help="MS1 정밀질량 게이트 ppm (정량 adduct 인정 기준, 기본 5)")
    p.add_argument("--max-charge", type=int, default=3, help="최대 전하수 (기본 3)")
    p.add_argument("--no-ms2", action="store_true",
                   help="MS2 근거 없이도 정량(기본은 MS2 확인된 글리칸만)")
    p.add_argument("--min-intensity", type=float, default=0.0,
                   help="이 강도 미만 adduct 무시 (기본 0)")
    p.add_argument("--keep-mzml", action="store_true", help="변환된 mzML 보존")
    # 조성 탐색 범위(선택)
    p.add_argument("--hexnac", nargs=2, type=int, metavar=("MIN", "MAX"))
    p.add_argument("--hex", nargs=2, type=int, metavar=("MIN", "MAX"))
    p.add_argument("--fuc", nargs=2, type=int, metavar=("MIN", "MAX"))
    p.add_argument("--neu5ac", nargs=2, type=int, metavar=("MIN", "MAX"))
    p.add_argument("--neu5gc", nargs=2, type=int, metavar=("MIN", "MAX"))
    p.add_argument("--xyl", action="store_true", help="Xyl 포함 탐색")
    args = p.parse_args(argv)

    if not os.path.isfile(args.input):
        print(f"[오류] 입력 파일이 없습니다: {args.input}")
        return 2

    t0 = time.time()
    out = args.output or (os.path.splitext(args.input)[0] + "_glycans.xlsx")
    work = os.path.join(os.path.dirname(os.path.abspath(args.input)), "_glycan_work")

    # 1) .raw -> mzML
    mzml_path = raw.to_mzml(args.input, work)

    # 2) 파싱
    data = mzml_parse.parse(mzml_path)

    # 3) 후보 생성
    ranges = build_ranges(args)
    cands = compositions.generate(ranges)
    print(f"[조성] 후보 {len(cands)}개 생성")

    # 4) 동정 + 정량
    results = identify.run(
        cands, msdata=data, max_charge=args.max_charge, ppm_tol=args.ppm,
        ms1_ppm=args.ms1_ppm, require_ms2=not args.no_ms2,
        min_intensity=args.min_intensity)

    if not results:
        print("[결과] 동정된 글리칸이 없습니다. --ppm 를 늘리거나 --no-ms2 를 시도하세요.")
        return 1

    # 5) 리포트
    meta = {"sample": os.path.basename(args.input), "source": "ThermoRawFileParser + glycan_tool"}
    report.write(results, out, meta=meta)

    # 정리
    if not args.keep_mzml and mzml_path.startswith(work):
        try:
            os.remove(mzml_path)
            os.rmdir(work)
        except OSError:
            pass

    dt = time.time() - t0
    print(f"\n[완료] 글리칸 {len(results)}개 -> {out}  ({dt:.1f}s)")
    print("상위 5개:")
    for r in results[:5]:
        print(f"  {r['name']:32s} {r['type']:12s} {r['relative_pct']:5.2f}%  "
              f"(m/z {r['best_mz']:.4f} {r['best_adduct']} {r['best_z']}+, MS2 {r['ms2_count']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
