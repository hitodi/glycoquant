"""
질량 엔진 (mass engine)
-----------------------
ProA로 환원말단 표지(reductive amination)된 N-글리칸의 monoisotopic
중성질량과 adduct 이온의 m/z를 계산한다.

검증: 기존 'Alzhemer of brain.xlsx > 계산기' 시트의 AC열(=[M+H]+) 44개를
     이 모델로 재계산한 결과 worst 0.0007 ppm 으로 일치함.

모델 정의 (계산기 시트 수식과 대수적으로 동일):
  - 중성질량 M = Σ(residue_i) + H2O + (ProA_free - O)
      residue_i = (자유 단당질량) - H2O
      ProA 환원아민화는 자유 ProA에서 산소 1개(O)를 제거한 것으로 처리
  - 중성질량은 '표준' monoisotopic 원소질량을 사용
  - 양이온(adduct)은 전자질량을 뺀 'bare cation' 질량을 사용 (H+, Na+, K+)
  - m/z = (M + Σ cation_mass) / (cation 개수=전하수 z)
"""

# 표준 monoisotopic 원소질량 (중성 원자)
ELEMENT = {
    "C": 12.0,
    "H": 1.0078250319,
    "N": 14.0030740052,
    "O": 15.9949146221,
    "S": 31.97207069,
    "Na": 22.98976928,
    "K": 38.96370649,
}

ELECTRON = 0.00054857991

# 자유(free) 단당 / 표지의 원소 조성 (계산기 시트 기준)
FREE_FORMULA = {
    "HexNAc": {"C": 8, "H": 15, "N": 1, "O": 6},
    "Hex":    {"C": 6, "H": 12, "O": 6},
    "dHex":   {"C": 6, "H": 12, "O": 5},   # = Fuc
    "Neu5Ac": {"C": 11, "H": 19, "N": 1, "O": 9},
    "Neu5Gc": {"C": 11, "H": 19, "N": 1, "O": 10},
    "Xyl":    {"C": 5, "H": 10, "O": 5},
}
PROA_FORMULA = {"C": 13, "H": 21, "N": 3, "O": 1}

# 단당 종류 (ProA 라벨은 별도 취급)
MONOSACCHARIDES = ("HexNAc", "Hex", "dHex", "Neu5Ac", "Neu5Gc", "Xyl")

# adduct 양이온 (bare cation = 중성원자 - 전자)
CATION_MASS = {
    "H": ELEMENT["H"] - ELECTRON,    # 1.0072765
    "Na": ELEMENT["Na"] - ELECTRON,  # 22.9892207
    "K": ELEMENT["K"] - ELECTRON,    # 38.9631579
}


def formula_mass(formula: dict) -> float:
    """원소 조성 dict -> monoisotopic 질량."""
    return sum(ELEMENT[el] * n for el, n in formula.items())


_H2O = formula_mass({"H": 2, "O": 1})
_O = ELEMENT["O"]


def neutral_mass(composition: dict) -> float:
    """
    글리칸 조성 -> ProA 표지된 중성 monoisotopic 질량.

    composition 예: {"HexNAc":2,"Hex":3,"dHex":1,"Neu5Ac":0,"Neu5Gc":0,"Xyl":0,"ProA":1}
    (없는 키는 0으로 간주)
    """
    res = 0.0
    for s in MONOSACCHARIDES:
        n = composition.get(s, 0)
        if n:
            res += (formula_mass(FREE_FORMULA[s]) - _H2O) * n
    proa = composition.get("ProA", 1)  # 라벨은 보통 1개
    return res + _H2O + (formula_mass(PROA_FORMULA) - _O) * proa


def element_formula(composition: dict) -> dict:
    """글리칸 조성 -> ProA 표지 후 최종 원소 조성(C/H/N/O)."""
    tot = {"C": 0, "H": 0, "N": 0, "O": 0}
    n_sug = 0
    for s in MONOSACCHARIDES:
        n = composition.get(s, 0)
        n_sug += n
        for el, c in FREE_FORMULA[s].items():
            tot[el] += c * n
    proa = composition.get("ProA", 1)
    for el, c in PROA_FORMULA.items():
        tot[el] += c * proa
    bonds = max(n_sug - 1, 0)            # 글리코사이드 결합 = (단당수 - 1)
    tot["H"] -= 2 * bonds                # 결합당 H2O 제거
    tot["O"] -= 1 * bonds
    tot["O"] -= 1 * proa                 # ProA 환원아민화 -O
    return {k: v for k, v in tot.items() if v}


def formula_str(composition: dict) -> str:
    f = element_formula(composition)
    return "".join(f"{el}{f[el]}" for el in ("C", "H", "N", "O", "S") if f.get(el))


def ion_mz(neutral: float, adduct, z: int) -> float:
    """
    중성질량 + adduct 조합 -> m/z.
    adduct: 양이온 종류의 시퀀스. 예: ("H",)=[M+H]+, ("H","Na")=[M+H+Na]2+
            len(adduct)가 곧 전하수 z 여야 한다.
    """
    if len(adduct) != z:
        raise ValueError("adduct 개수와 전하수 z가 일치해야 합니다.")
    return (neutral + sum(CATION_MASS[a] for a in adduct)) / z
