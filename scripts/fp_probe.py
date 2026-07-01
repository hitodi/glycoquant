"""
알츠하이머 FULL-후보 위양성(FP) 프로브
--------------------------------------
알츠하이머 시료는 시알산이 거의 0% (보고용 시트 확인) 이므로,
FULL 조성후보 + MS1-first 로 돌렸을 때 '컷오프를 넘는 시알산 글리칸'은
전부 정의상 위양성이다. 그 개수가 곧 MS1-first 의 정밀도 지표.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco import config as cfgmod, compositions, mzml_parse, identify
from glyco.chem import Chemistry

MZML = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "tests", "data", "sample.mzML")
_D = {}


def probe(ms1_first, noise_factor=15.0, min_adducts=3, ppm=3.0, out_cutoff=0.5):
    cfg = cfgmod.load(); chem = Chemistry(cfg)
    cands = compositions.from_config(cfg, chem=chem)   # FULL 조합 후보
    if "d" not in _D:
        _D["d"] = mzml_parse.parse(MZML, keep_ms2_peaks=True, log=lambda *a: None)
    data = _D["d"]
    res = identify.run(cands, msdata=data, ppm_tol=10, ms1_ppm=5, require_ms2=True,
                       quant_method="area", rt_consistency=1.0, chem=chem,
                       require_diagnostic=cfg.require_diagnostic, ms2_ppm=cfg.ms2_ppm,
                       ms1_first=ms1_first, ms1_min_adducts=min_adducts,
                       ms1_first_ppm=ppm, ms1_noise_factor=noise_factor, log=lambda *a: None)
    tot = sum(r["intensity_sum"] for r in res) or 1
    for r in res:
        r["pct"] = r["intensity_sum"] / tot * 100
    over = [r for r in res if r["pct"] >= out_cutoff]
    sia = [r for r in over if r["composition"].get("Neu5Ac", 0) or r["composition"].get("Neu5Gc", 0)]
    lab = f"ms1_first={ms1_first}" + (f" noise×{noise_factor:g} ≥{min_adducts} {ppm}ppm" if ms1_first else "")
    print(f"[{lab}] 동정 {len(res)} | ≥{out_cutoff}% {len(over)}개 | "
          f"그중 시알산(=위양성) {len(sia)}개: {[r['oxford'] for r in sia][:8]}")
    return len(sia)


if __name__ == "__main__":
    probe(False)
    for rf in (1e-4, 3e-4, 1e-3):
        probe(True, rel_floor=rf, min_adducts=3, ppm=3.0)
