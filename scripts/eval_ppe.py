"""
논문 PPE(Table 1) 정답 대비 평가 하네스
---------------------------------------
pancreatin 시료를 도구로 분석 -> 조성으로 논문 32개와 매칭 ->
recall / Pearson / Spearman / 중앙오차 / NeuGc% / FA2% / 유형분포 를 출력.
식(파라미터)을 바꿔가며 이 지표가 오르는지 보는 용도.

사용:  python scripts/eval_ppe.py            # 기본(출하) 설정
       evaluate(overrides=...) 로 파라미터 실험
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco import config as cfgmod, compositions, mzml_parse, identify
from glyco.chem import Chemistry

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MZML = os.path.join(ROOT, "tests", "data", "pancreatin.mzML")
TRUTH = os.path.join(ROOT, "tests", "data", "paper_ppe.json")

_CACHE = {}


def _key(c):
    return (c["HexNAc"], c["Hex"], c["dHex"], c.get("Neu5Ac", 0), c.get("Neu5Gc", 0))


def _pearson(a, b):
    n = len(a); ma = sum(a) / n; mb = sum(b) / n
    va = math.sqrt(sum((x - ma) ** 2 for x in a)); vb = math.sqrt(sum((x - mb) ** 2 for x in b))
    return sum((a[i] - ma) * (b[i] - mb) for i in range(n)) / (va * vb) if va * vb else 0.0


def _spearman(a, b):
    n = len(a)
    def rk(x):
        o = sorted(range(n), key=lambda i: x[i]); r = [0] * n
        for p, i in enumerate(o): r[i] = p
        return r
    ra, rb = rk(a), rk(b)
    return 1 - 6 * sum((ra[i] - rb[i]) ** 2 for i in range(n)) / (n * (n * n - 1))


def _load_data():
    if "data" not in _CACHE:
        _CACHE["data"] = mzml_parse.parse(MZML, keep_ms2_peaks=True, log=lambda *a: None)
    return _CACHE["data"]


def evaluate(config_path=None, overrides=None, label="", verbose=True):
    overrides = overrides or {}
    cfg = cfgmod.load(config_path)
    chem = Chemistry(cfg)
    truth = json.load(open(TRUTH))
    truth_by = {_key(t["composition"]): t for t in truth}

    cands = compositions.from_config(cfg, chem=chem)
    data = _load_data()
    q = cfg.quantify
    res = identify.run(
        cands, msdata=data,
        ppm_tol=overrides.get("precursor_ppm", cfg.precursor_ppm),
        ms1_ppm=overrides.get("ms1_ppm", cfg.ms1_ppm),
        max_charge=overrides.get("max_charge", cfg.max_charge),
        require_ms2=overrides.get("require_ms2", q.get("require_ms2", True)),
        quant_method=overrides.get("quant_method", q.get("method", "area")),
        rt_window=overrides.get("rt_window", q.get("rt_window_min", 0.5)),
        rt_consistency=overrides.get("rt_consistency", q.get("rt_consistency_min", 1.0)),
        min_intensity=overrides.get("min_intensity", 0.0),
        chem=chem,
        require_diagnostic=(None if overrides.get("no_diagnostic")
                            else overrides.get("require_diagnostic", cfg.require_diagnostic)),
        ms2_ppm=overrides.get("ms2_ppm", cfg.ms2_ppm),
        ms1_first=overrides.get("ms1_first", False),
        ms1_min_adducts=overrides.get("ms1_min_adducts", 3),
        ms1_first_ppm=overrides.get("ms1_first_ppm", 5.0),
        ms1_noise_factor=overrides.get("ms1_noise_factor", 15.0),
        log=lambda *a: None,
    )
    tot = sum(r["intensity_sum"] for r in res) or 1.0
    for r in res:
        r["pct"] = r["intensity_sum"] / tot * 100
    res_by = {_key(r["composition"]): r for r in res}

    # 논문과 같은 기준(≥0.6% 컷오프 후 재정규화)으로 헤드라인 지표 계산
    cut = [r for r in res if r["pct"] >= 0.6]
    ctot = sum(r["intensity_sum"] for r in cut) or 1.0
    cut_by = {}
    for r in cut:
        r["cpct"] = r["intensity_sum"] / ctot * 100
        cut_by[_key(r["composition"])] = r
    from collections import Counter
    ctype = Counter()
    for r in cut:
        ctype[r["type"]] += r["cpct"]
    neugc_cut = sum(r["cpct"] for r in cut if r["composition"].get("Neu5Gc", 0))
    fa2_cut = cut_by.get((4, 3, 1, 0, 0), {}).get("cpct", 0)
    fa2g1sg1 = cut_by.get((4, 4, 1, 0, 1), {}).get("cpct", 0)  # 논문 #4 = 6.6%

    # 매칭 (논문 조성 기준)
    shared = [k for k in truth_by if k in res_by]
    tp = [truth_by[k]["pct"] for k in shared]
    op = [res_by[k]["pct"] for k in shared]
    # 공유집합 재정규화(분모 일치)
    st, so = sum(tp), sum(op)
    tpn = [x / st * 100 for x in tp]; opn = [x / so * 100 for x in op]
    aerr = sorted(abs(tpn[i] - opn[i]) for i in range(len(shared)))
    metrics = {
        "n_detected": len(res),
        "recall": f"{len(shared)}/{len(truth)}",
        "recall_n": len(shared),
        "pearson": _pearson(tpn, opn) if len(shared) > 2 else 0,
        "spearman": _spearman(tpn, opn) if len(shared) > 2 else 0,
        "med_err": aerr[len(aerr) // 2] if aerr else 0,
        "max_err": max(aerr) if aerr else 0,
        # 관심 지표
        "NeuGc_pct": sum(r["pct"] for r in res if r["composition"].get("Neu5Gc", 0)),
        "FA2_pct": res_by.get((4, 3, 1, 0, 0), {}).get("pct", 0),
        "missing": [truth_by[k]["name"] for k in truth_by if k not in res_by],
    }
    metrics.update({"n_cut": len(cut), "NeuGc_cut": neugc_cut, "FA2_cut": fa2_cut,
                    "FA2G1Sg1_cut": fa2g1sg1,
                    "type_cut": {t: round(ctype[t], 1) for t in ("High-mannose", "Hybrid", "Complex")}})
    if verbose:
        print(f"\n[{label or 'eval'}] 동정 {metrics['n_detected']}개 (≥0.6% {len(cut)}개, 논문 32)"
              f" | recall {metrics['recall']} | Spearman {metrics['spearman']:+.3f}")
        print(f"        [≥0.6% 기준] NeuGc {neugc_cut:.1f}%(논문14.0) | FA2 {fa2_cut:.1f}%(논문17.7)"
              f" | FA2G1Sg1 {fa2g1sg1:.1f}%(논문6.6)")
        print(f"        유형 {metrics['type_cut']} (논문 HM31.3/Hyb4.3/Cplx64.4)"
              f" | 놓침 {len(metrics['missing'])}: {', '.join(metrics['missing'][:10])}")
    return metrics


if __name__ == "__main__":
    evaluate(label="baseline (출하 기본)")
