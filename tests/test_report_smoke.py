"""리포트 writer 스모크: write / write_diagnostic_screening / write_aggregated 가
_style_body 등 공용 헬퍼를 거쳐 xlsx 를 예외 없이 만들어내는지 확인(내용 검증 아님).
이 writer 들은 기존 스위트에 커버가 없어, 스타일 헬퍼 리팩토링의 안전망으로 추가."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from glyco import report, targets


def _result(oxford, comp, typ, intensity):
    return {
        "oxford": oxford, "name": oxford, "type": typ,
        "composition": comp,
        "sialylated": comp.get("Neu5Ac", 0) > 0, "fucosylated": comp.get("dHex", 0) > 0,
        "best_mz": 800.4, "best_adduct": "H", "best_z": 2, "rt": 12.3,
        "ms2_count": 3, "evidence": "MS2", "best_ppm": 1.2, "intensity_sum": intensity,
        "relative_pct": 0.0,
        "adducts": [{"adduct": "H", "z": 2, "mz": 800.4, "rt": 12.3, "intensity": intensity}],
    }


def _results():
    return [
        _result("M5", {"HexNAc": 2, "Hex": 5, "dHex": 0, "Neu5Ac": 0, "Neu5Gc": 0}, "High-mannose", 5000.0),
        _result("FA2", {"HexNAc": 4, "Hex": 3, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 0}, "Complex", 3000.0),
    ]


def test_write_glycans(tmp_path):
    out = str(tmp_path / "g.xlsx")
    report.write(_results(), out, meta={"sample": "s.raw", "source": "test"})
    assert os.path.exists(out)


def test_write_aggregated(tmp_path):
    from glyco import aggregate
    agg = aggregate.aggregate([("rep1", _results()), ("rep2", _results())])
    out = str(tmp_path / "a.xlsx")
    report.write_aggregated(agg, out, meta={"group": "grp"})
    assert os.path.exists(out)


def test_write_diagnostic_screening(tmp_path):
    from glyco.identify import screen_diagnostics
    spec = targets.DiagnosticSpec(
        ions={"ProA-HexNAc": 441.2708, "HexNAc": 204.0867},
        accept=[{"any": ["ProA-HexNAc"], "min": 1}, {"any": ["HexNAc"], "min": 1}],
        features=[{"name": "Neu5Ac", "mz": [292.1027]}],
        ppm=5.0, group_ppm=5.0, split_by_charge=True)

    class _FakeMs:
        ms2 = [(1., 500., 2, "A"), (2., 500.001, 2, "B")]
        ms2_peaks = {"A": (np.array([441.2708, 204.0867, 292.1027]), np.array([9., 9., 9.])),
                     "B": (np.array([441.2708, 204.0867]), np.array([9., 9.]))}
        ms2_info = {}
    rows = screen_diagnostics(_FakeMs(), spec)
    out = str(tmp_path / "d.xlsx")
    report.write_diagnostic_screening(rows, spec, out, meta={"sample": "s.raw"})
    assert os.path.exists(out) and len(rows) == 2


if __name__ == "__main__":
    import pathlib
    import tempfile
    d = pathlib.Path(tempfile.mkdtemp())
    test_write_glycans(d)
    test_write_aggregated(d)
    test_write_diagnostic_screening(d)
    print("report 스모크 통과 ✓")
