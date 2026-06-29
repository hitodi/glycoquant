"""Oxford 명명법이 논문 Table 1 예시와 일치하는지."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco.nomenclature import oxford

CASES = [
    ({"HexNAc": 2, "Hex": 5, "dHex": 0, "Neu5Ac": 0, "Neu5Gc": 0}, "M5"),
    ({"HexNAc": 2, "Hex": 3, "dHex": 0, "Neu5Ac": 0, "Neu5Gc": 0}, "M3"),
    ({"HexNAc": 2, "Hex": 3, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 0}, "FM3"),
    ({"HexNAc": 4, "Hex": 3, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 0}, "FA2"),
    ({"HexNAc": 4, "Hex": 5, "dHex": 1, "Neu5Ac": 2, "Neu5Gc": 0}, "FA2G2S2"),
    ({"HexNAc": 4, "Hex": 4, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 1}, "FA2G1Sg1"),
    ({"HexNAc": 4, "Hex": 3, "dHex": 0, "Neu5Ac": 0, "Neu5Gc": 0}, "A2"),
    ({"HexNAc": 3, "Hex": 3, "dHex": 1, "Neu5Ac": 0, "Neu5Gc": 0}, "FA1"),
]


def test_oxford_matches_paper():
    for comp, expected in CASES:
        assert oxford(comp) == expected, f"{comp} -> {oxford(comp)} (기대 {expected})"


if __name__ == "__main__":
    test_oxford_matches_paper()
    print("Oxford 명명법 테스트 통과 ✓")
