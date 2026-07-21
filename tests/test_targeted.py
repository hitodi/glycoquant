"""타깃 스크리닝: floor 그룹핑 경계 + K-of-N 티어 + 타깃파일 검증 (mzML 불필요)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest

from glyco.identify import floor_bin, targeted_screen


def test_floor_bin_boundary():
    # 같은 0.1 구간은 한 그룹, 경계 넘으면 갈림
    assert abs(floor_bin(557.70) - 557.7) < 1e-9
    assert abs(floor_bin(557.74) - 557.7) < 1e-9
    assert abs(floor_bin(557.79) - 557.7) < 1e-9
    assert abs(floor_bin(557.69) - 557.6) < 1e-9
    assert abs(floor_bin(557.80) - 557.8) < 1e-9


class _FakeMs:
    def __init__(self, ms2, peaks):
        self.ms2 = ms2                 # [(rt,pmz,z,scan)]
        self.ms2_peaks = peaks         # {scan:(mz,inten)}


def test_k_of_n_tiering():
    glycans = [{"name": "G", "ions": [100.0, 200.0, 300.0]}]   # N=3
    peaks = {
        "1": (np.array([100.0, 200.0]), np.array([9., 9.])),   # 2 hit → matched
        "2": (np.array([100.0]), np.array([9.])),              # 1 hit → holding
        "3": (np.array([999.0]), np.array([9.])),              # 0 hit → 제외
    }
    ms2 = [(10.0, 500.0, 2, "1"), (11.0, 501.0, 2, "2"), (12.0, 502.0, 2, "3")]
    rows = targeted_screen(_FakeMs(ms2, peaks), glycans, ppm=20, min_hits=2)
    by = {r["scan"]: r for r in rows}
    assert by["1"]["tier"] == "matched" and by["1"]["n_hit"] == 2
    assert by["2"]["tier"] == "holding" and by["2"]["n_hit"] == 1
    assert "3" not in by                                       # 0 hit → 행 없음


def test_multi_glycan_per_pair_tier():
    """한 스캔이 글리칸별로 다른 티어 가능(per-pair 판정)."""
    glycans = [{"name": "A", "ions": [100.0, 200.0, 300.0]},
               {"name": "B", "ions": [100.0, 700.0, 800.0]}]
    peaks = {"1": (np.array([100.0, 200.0, 300.0]), np.array([9., 9., 9.]))}  # A=3, B=1
    rows = targeted_screen(_FakeMs([(10.0, 500.0, 2, "1")], peaks), glycans, ppm=20, min_hits=2)
    by = {r["glycan"]: r for r in rows}
    assert by["A"]["tier"] == "matched"      # 3/3
    assert by["B"]["tier"] == "holding"      # 1/3 (공유 이온 100만)


def test_targets_validation(tmp_path):
    from glyco import targets
    good = tmp_path / "t.yaml"
    good.write_text("min_hits: 2\nglycans:\n  - name: A\n    ions: [204.0867, 366.1395, 528.19]\n", encoding="utf-8")
    spec = targets.load(str(good))
    assert len(spec.glycans) == 1 and spec.min_hits == 2
    # 이온 수 < min_hits → 거부
    bad = tmp_path / "b.yaml"
    bad.write_text("min_hits: 3\nglycans:\n  - name: A\n    ions: [204.0867]\n", encoding="utf-8")
    with pytest.raises(ValueError):
        targets.load(str(bad))
    # 숫자 아닌 이온 → 거부
    bad2 = tmp_path / "b2.yaml"
    bad2.write_text("glycans:\n  - name: A\n    ions: [abc, 204.0]\n", encoding="utf-8")
    with pytest.raises(ValueError):
        targets.load(str(bad2))


if __name__ == "__main__":
    test_floor_bin_boundary(); test_k_of_n_tiering(); test_multi_glycan_per_pair_tier()
    print("타깃 스크리닝 테스트 통과 ✓")
