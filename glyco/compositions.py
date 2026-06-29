"""
조성 생성기 + adduct 이온 목록 (설정 구동)
------------------------------------------
새 .raw 에는 기존과 다른 글리칸이 있을 수 있으므로, 설정의 탐색범위에서
조성을 '조합으로 생성'한다(하드코딩 금지). 단당 종류·범위·타당성 규칙은 Config 에서 온다.
"""

from itertools import combinations_with_replacement, product
from dataclasses import dataclass, field

from . import masses
from .nomenclature import oxford


# --- 하위호환용 기본 범위(설정 없이 호출 시) ---
@dataclass
class SearchRanges:
    HexNAc: range = field(default_factory=lambda: range(2, 9))
    Hex: range = field(default_factory=lambda: range(3, 13))
    dHex: range = field(default_factory=lambda: range(0, 5))
    Neu5Ac: range = field(default_factory=lambda: range(0, 5))
    Neu5Gc: range = field(default_factory=lambda: range(0, 5))
    Xyl: range = field(default_factory=lambda: range(0, 1))
    ProA: int = 1

    def as_dict(self):
        return {k: getattr(self, k) for k in
                ("HexNAc", "Hex", "dHex", "Neu5Ac", "Neu5Gc", "Xyl")}


def _plausible(c: dict, rules: dict) -> bool:
    """N-글리칸 타당성 휴리스틱 (설정 토글)."""
    hexnac = c.get("HexNAc", 0)
    hexn = c.get("Hex", 0)
    sia = c.get("Neu5Ac", 0) + c.get("Neu5Gc", 0)
    fuc = c.get("dHex", 0)
    if hexnac < rules.get("min_hexnac", 2):
        return False
    if hexn < rules.get("min_hex", 3):
        return False
    antennae = max(hexnac - 2, 0)
    if rules.get("sialic_le_antennae", True) and sia > antennae:
        return False
    if rules.get("fucose_le_antennae_plus1", True) and fuc > antennae + 1:
        return False
    if rules.get("highmannose_no_sia_fuc", True) and hexnac == 2 and (sia or fuc):
        return False
    return True


def generate(ranges=None, plausible_only=True, rules=None, label_count=1):
    """
    조성 생성. ranges 는 {name: range} dict 또는 SearchRanges.
    각 결과: {composition, neutral, formula, name, oxford}
    """
    if ranges is None:
        ranges = SearchRanges()
    rdict = ranges.as_dict() if isinstance(ranges, SearchRanges) else dict(ranges)
    rules = rules or {}
    names = list(rdict.keys())
    out = []
    for combo in product(*[list(rdict[n]) for n in names]):
        c = dict(zip(names, combo))
        c["ProA"] = label_count
        if plausible_only and not _plausible(c, rules):
            continue
        out.append({
            "composition": c,
            "neutral": masses.neutral_mass(c),
            "formula": masses.formula_str(c),
            "name": composition_name(c),
            "oxford": oxford(c),
        })
    return out


def from_config(cfg, plausible_only=True):
    """Config 의 search_ranges/plausibility 로 조성 생성."""
    ranges = {n: range(lo, hi + 1) for n, (lo, hi) in cfg.search_ranges.items()}
    return generate(ranges=ranges, plausible_only=plausible_only,
                    rules=cfg.plausibility, label_count=cfg.label.count)


def composition_name(c: dict) -> str:
    parts = []
    for key, label in (("HexNAc", "HexNAc"), ("Hex", "Hex"), ("dHex", "Fuc"),
                       ("Neu5Ac", "Neu5Ac"), ("Neu5Gc", "Neu5Gc"), ("Xyl", "Xyl")):
        if c.get(key, 0):
            parts.append(f"{label}{c[key]}")
    return " ".join(parts) if parts else "(empty)"


def adduct_ions(max_charge=3, cations=("H", "Na", "K")):
    """전하 1..max_charge 의 모든 adduct 다중집합. [(adduct_tuple, z), ...]"""
    out = []
    for z in range(1, max_charge + 1):
        for combo in combinations_with_replacement(cations, z):
            out.append((combo, z))
    return out


def adduct_label(adduct) -> str:
    from collections import Counter
    cnt = Counter(adduct)
    parts = []
    for cat in ("H", "Na", "K"):
        if cnt[cat]:
            parts.append(f"{cnt[cat] if cnt[cat] > 1 else ''}{cat}")
    return "+".join(parts)
