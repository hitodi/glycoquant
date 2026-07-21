"""진단 스크리닝: 채택규칙(441 필수+204/366) · feature OR · 그룹핑 · 규칙파일 검증 (mzML 불필요)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest

from glyco.identify import screen_diagnostics, floor_bin
from glyco.report import group_by_precursor
from glyco import targets


class _FakeMs:
    def __init__(self, ms2, peaks, info=None):
        self.ms2 = ms2
        self.ms2_peaks = peaks
        self.ms2_info = info or {}


def _spec():
    return targets.DiagnosticSpec(
        ions={"ProA-HexNAc": 441.2708, "HexNAc": 204.0867, "HexNAc+Hex": 366.1395},
        accept=[{"any": ["ProA-HexNAc"], "min": 1},
                {"any": ["HexNAc", "HexNAc+Hex"], "min": 1}],
        features=[{"name": "Neu5Ac", "mz": [292.1027]},
                  {"name": "bisecting", "mz": [1009.4823, 1155.5403]}],
        ppm=5.0, group_ppm=5.0, split_by_charge=True)


def test_accept_rule_441_required():
    spec = _spec()
    peaks = {
        "A": (np.array([441.2708, 204.0867]), np.array([9., 9.])),   # 441+204 → 채택
        "B": (np.array([441.2708, 366.1395]), np.array([9., 9.])),   # 441+366 → 채택
        "C": (np.array([204.0867, 366.1395]), np.array([9., 9.])),   # 441 없음 → 버림
        "D": (np.array([441.2708]), np.array([9.])),                 # 441만(204/366 없음) → 버림
    }
    ms2 = [(1., 500., 2, "A"), (2., 501., 2, "B"), (3., 502., 2, "C"), (4., 503., 2, "D")]
    rows = screen_diagnostics(_FakeMs(ms2, peaks), spec)
    got = {r["scan"] for r in rows}
    assert got == {"A", "B"}          # 441 + (204 or 366) 만 채택


def test_feature_annotation_or():
    spec = _spec()
    peaks = {"A": (np.array([441.2708, 204.0867, 292.1027, 1155.5403]), np.array([9.] * 4))}
    rows = screen_diagnostics(_FakeMs([(1., 500., 2, "A")], peaks), spec)
    feats = set(rows[0]["features"])
    assert "Neu5Ac" in feats               # 292 검출
    assert "bisecting" in feats            # 1155(OR) 검출


def test_monoisotope_from_info():
    spec = _spec()
    peaks = {"A": (np.array([441.2708, 204.0867]), np.array([9., 9.]))}
    info = {"A": {"monoisotopic_mz": 765.3788}}
    rows = screen_diagnostics(_FakeMs([(1., 766.38, 1, "A")], peaks, info), spec)
    assert abs(rows[0]["monoisotope"] - 765.3788) < 1e-6   # userParam 사용
    # info 없으면 precursor 폴백
    rows2 = screen_diagnostics(_FakeMs([(1., 766.38, 1, "A")], peaks), spec)
    assert abs(rows2[0]["monoisotope"] - 766.38) < 1e-6


def test_grouping_precursor_and_charge():
    rows = [{"precursor": 500.0000, "charge": 2, "scan": "1"},
            {"precursor": 500.0015, "charge": 2, "scan": "2"},   # 3ppm → 같은 그룹
            {"precursor": 500.0100, "charge": 2, "scan": "3"},   # 20ppm → 다른 그룹
            {"precursor": 500.0000, "charge": 3, "scan": "4"}]   # charge 다름 → 다른 그룹
    g = dict((r["scan"], gid) for gid, r in group_by_precursor(rows, ppm=5, split_by_charge=True))
    assert g["1"] == g["2"]
    assert g["3"] != g["1"]
    assert g["4"] != g["1"]


def test_floor_bin_boundary():
    assert abs(floor_bin(557.74) - 557.7) < 1e-9
    assert abs(floor_bin(557.69) - 557.6) < 1e-9


def test_rule_validation(tmp_path):
    good = tmp_path / "d.yaml"
    good.write_text(
        "ppm: 5\nions:\n  A: 441.2708\n  B: 204.0867\n"
        "accept:\n  - {any: [A], min: 1}\n  - {any: [B], min: 1}\n", encoding="utf-8")
    spec = targets.load(str(good))
    assert spec.ppm == 5 and len(spec.ions) == 2
    # accept 가 ions 에 없는 이름 참조 → 거부
    bad = tmp_path / "b.yaml"
    bad.write_text("ions:\n  A: 441.0\naccept:\n  - {any: [ZZZ], min: 1}\n", encoding="utf-8")
    with pytest.raises(ValueError):
        targets.load(str(bad))


if __name__ == "__main__":
    for fn in [test_accept_rule_441_required, test_feature_annotation_or,
               test_monoisotope_from_info, test_grouping_precursor_and_charge, test_floor_bin_boundary]:
        fn()
    print("진단 스크리닝 테스트 통과 ✓")
