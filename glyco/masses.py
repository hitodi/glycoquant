"""
질량 엔진 (하위호환 shim)
-------------------------
실제 계산은 chem.Chemistry(설정 구동)가 한다. 이 모듈은 기존 코드/테스트가
쓰던 모듈수준 API(neutral_mass/ion_mz/formula_mass/FREE_FORMULA/CATION_MASS)를
'기본 ProA 설정'으로 그대로 제공한다.

검증: 계산기 시트 44개 [M+H]+ < 0.001 ppm (tests/test_masses.py).
모델: 중성 M = Σ(residue) + H2O + (label_free - attach_loss);
      표준 monoisotopic 원소질량, cation 은 전자 뺀 bare 질량.
"""

from . import config as _config
from .chem import Chemistry, ELEMENT, ELECTRON, formula_mass

# 기본(ProA) 설정으로 만든 전역 화학 엔진
_default = Chemistry(_config.load())

# --- 하위호환 상수 ---
FREE_FORMULA = _default.monos          # 단당 자유조성 dict
PROA_FORMULA = _default.label.formula
MONOSACCHARIDES = _default.mono_names
CATION_MASS = _default.cation          # {'H':..,'Na':..,'K':..}


def neutral_mass(composition: dict) -> float:
    return _default.neutral_mass(composition)


def element_formula(composition: dict) -> dict:
    return _default.element_formula(composition)


def formula_str(composition: dict) -> str:
    return _default.formula_str(composition)


def ion_mz(neutral: float, adduct, z: int) -> float:
    return _default.ion_mz(neutral, adduct, z)


__all__ = ["ELEMENT", "ELECTRON", "formula_mass", "FREE_FORMULA", "PROA_FORMULA",
           "MONOSACCHARIDES", "CATION_MASS", "neutral_mass", "element_formula",
           "formula_str", "ion_mz", "Chemistry"]
