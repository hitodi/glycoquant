"""반복 취합(aggregate) 로직 단위 테스트 — mzML 불필요(합성 데이터)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco.aggregate import aggregate


def _g(hexnac, hexn, dhex, ac, gc, pct, typ="Complex", evidence="MS2"):
    return {"composition": {"HexNAc": hexnac, "Hex": hexn, "dHex": dhex,
                            "Neu5Ac": ac, "Neu5Gc": gc},
            "oxford": "X", "name": "x", "type": typ,
            "sialylated": bool(ac or gc), "fucosylated": bool(dhex),
            "relative_pct": pct, "evidence": evidence}


def test_identical_replicates_zero_sd():
    """동일 반복 → 평균=값, SD=0, n/N=2/2."""
    r = [_g(4, 3, 1, 0, 0, 20.0)]
    agg = aggregate([("rep1", r), ("rep2", r)])
    g = agg["glycans"][0]
    assert g["mean"] == 20.0 and g["sd"] == 0.0
    assert g["n_detected"] == 2 and g["n_total"] == 2


def test_mean_sd_and_partial_detection():
    """값이 다르면 평균/SD 계산, 한쪽만 검출되면 n/N=1/2."""
    a = [_g(4, 3, 1, 0, 0, 10.0), _g(2, 5, 0, 0, 0, 5.0)]   # FA2형 + M5형
    b = [_g(4, 3, 1, 0, 0, 20.0)]                            # FA2형만
    agg = aggregate([("A", a), ("B", b)])
    by = {tuple(sorted(g["composition"].items())): g for g in agg["glycans"]}
    fa2 = by[tuple(sorted({"HexNAc": 4, "Hex": 3, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 0}.items()))]
    m5 = by[tuple(sorted({"HexNAc": 2, "Hex": 5, "dHex": 0, "Neu5Ac": 0, "Neu5Gc": 0}.items()))]
    assert fa2["mean"] == 15.0 and fa2["n_detected"] == 2       # (10+20)/2
    assert abs(fa2["sd"] - 7.0710678) < 1e-4                    # 표본 SD
    # 부분 검출: 결측=0 정책 → 값 [5,0] → 평균 2.5, SD 3.54, n/N=1/2
    assert m5["n_detected"] == 1 and m5["n_total"] == 2
    assert m5["mean"] == 2.5                                    # (5+0)/2, 결측=0
    assert m5["sd"] > 0                                         # 산발검출은 SD로 드러남(false precision 방지)


def test_evidence_merge():
    """반복마다 evidence 다르면 병합 표기."""
    agg = aggregate([("A", [_g(4, 3, 1, 0, 0, 10, evidence="MS2")]),
                     ("B", [_g(4, 3, 1, 0, 0, 12, evidence="MS1")])])
    assert set(agg["glycans"][0]["evidence"].split("+")) == {"MS1", "MS2"}


if __name__ == "__main__":
    test_identical_replicates_zero_sd()
    test_mean_sd_and_partial_detection()
    test_evidence_merge()
    print("aggregate 테스트 통과 ✓")
