"""
질량 엔진 회귀 테스트 (리팩터 가드레일)
---------------------------------------
기존 'Alzhemer of brain.xlsx > 계산기' 시트의 AC열([M+H]+) 44개를
ground truth 로 삼아, 어떤 리팩터 후에도 < 0.001 ppm 으로 재현되는지 검증.

실행:  pytest tests/test_masses.py     또는    python tests/test_masses.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco import masses

FIXTURE = os.path.join(os.path.dirname(__file__), "data", "calculator_fixture.json")
TOL_PPM = 0.001


def _load():
    return [r for r in json.load(open(FIXTURE)) if r.get("AC")]


def test_neutral_and_protonated_mass():
    """44개 글리칸의 [M+H]+ 가 시트 AC값과 < 0.001 ppm 일치."""
    worst = 0.0
    for r in _load():
        comp = {k: r[k] for k in ("HexNAc", "Hex", "dHex", "Neu5Ac", "Neu5Gc", "Xyl", "ProA")}
        calc = masses.ion_mz(masses.neutral_mass(comp), ("H",), 1)
        ppm = abs(calc - r["AC"]) / r["AC"] * 1e6
        worst = max(worst, ppm)
        assert ppm < TOL_PPM, f"glycan #{r['no']} ppm={ppm:.4f} (calc {calc} vs AC {r['AC']})"
    print(f"[OK] {len(_load())}개 검증, worst |ppm| = {worst:.5f}")


def test_diagnostic_oxonium_ions():
    """논문(§2.4)이 명시한 진단 oxonium/단편 이온 m/z 재현 (<2 ppm)."""
    # (이름, 이론 m/z, 단편 조성, ProA 포함여부)
    cases = [
        ("HexNAc",            204.0867, {"HexNAc": 1}, False),
        ("HexNAc+Hex",        366.1395, {"HexNAc": 1, "Hex": 1}, False),
        ("Neu5Ac",            292.1027, {"Neu5Ac": 1}, False),
        ("Neu5Gc",            308.0976, {"Neu5Gc": 1}, False),
        ("ProA-HexNAc",       441.2708, {"HexNAc": 1}, True),
        ("ProA-Fuc-HexNAc",   587.3281, {"HexNAc": 1, "dHex": 1}, True),
    ]
    for name, theo, frag, has_proa in cases:
        if has_proa:
            # 환원말단 단편: neutral_mass(라벨 포함) + proton
            comp = dict(frag); comp["ProA"] = 1
            mz = masses.ion_mz(masses.neutral_mass(comp), ("H",), 1)
        else:
            # oxonium(B형): 잔기질량 합 + proton (탈수형)
            res = sum((masses.formula_mass(masses.FREE_FORMULA[m]) -
                       masses.formula_mass({"H": 2, "O": 1})) * n
                      for m, n in frag.items())
            mz = res + masses.CATION_MASS["H"]
        ppm = abs(mz - theo) / theo * 1e6
        assert ppm < 2.0, f"{name}: calc {mz:.4f} vs theo {theo} ({ppm:.2f} ppm)"
        print(f"[OK] {name:18s} calc={mz:.4f} theo={theo} ({ppm:+.2f} ppm)")


if __name__ == "__main__":
    test_neutral_and_protonated_mass()
    test_diagnostic_oxonium_ions()
    print("\n모든 질량 테스트 통과 ✓")
