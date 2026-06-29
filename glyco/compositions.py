"""
조성 생성기 + adduct 이온 목록
------------------------------
새 .raw 파일에는 기존 45개와 다른 글리칸이 들어있을 수 있으므로,
N-글리칸으로 말이 되는 범위에서 조성을 '조합으로 생성'한다.
(45개 하드코딩 금지 — advisor 지침)
"""

from itertools import combinations_with_replacement
from dataclasses import dataclass, field

from . import masses


@dataclass
class SearchRanges:
    """조성 탐색 범위 (CLI/설정으로 조절 가능)."""
    HexNAc: range = field(default_factory=lambda: range(2, 9))   # 코어 2 ~ 8
    Hex: range = field(default_factory=lambda: range(3, 13))     # 코어 3 ~ high-mannose
    dHex: range = field(default_factory=lambda: range(0, 5))     # Fuc 0 ~ 4
    Neu5Ac: range = field(default_factory=lambda: range(0, 5))
    Neu5Gc: range = field(default_factory=lambda: range(0, 5))
    Xyl: range = field(default_factory=lambda: range(0, 1))      # 기본 끔
    ProA: int = 1


def _is_plausible_nglycan(c: dict) -> bool:
    """N-글리칸 조성 타당성 휴리스틱 (GlycoMod 류)."""
    hexnac, hexn = c["HexNAc"], c["Hex"]
    sia = c["Neu5Ac"] + c["Neu5Gc"]
    fuc = c["dHex"]
    # 최소 trimannosyl 코어: HexNAc>=2, Hex>=3 (또는 paucimannose 허용 시 완화)
    if hexnac < 2 or hexn < 3:
        return False
    # 안테나(비코어 HexNAc) 수
    antennae = hexnac - 2
    # 시알산은 안테나 수를 넘을 수 없음
    if sia > max(antennae, 0):
        return False
    # 푸코스는 (코어1 + 안테나당 1) 정도로 제한
    if fuc > antennae + 1:
        return False
    # high-mannose(HexNAc==2)는 시알/푸코 없음
    if hexnac == 2 and (sia > 0 or fuc > 0):
        return False
    return True


def generate(ranges: SearchRanges = None, plausible_only: bool = True):
    """타당한 글리칸 조성들을 생성 (각 조성에 neutral mass 포함)."""
    r = ranges or SearchRanges()
    out = []
    for hexnac in r.HexNAc:
        for hexn in r.Hex:
            for dhex in r.dHex:
                for ac in r.Neu5Ac:
                    for gc in r.Neu5Gc:
                        for xyl in r.Xyl:
                            c = {"HexNAc": hexnac, "Hex": hexn, "dHex": dhex,
                                 "Neu5Ac": ac, "Neu5Gc": gc, "Xyl": xyl,
                                 "ProA": r.ProA}
                            if plausible_only and not _is_plausible_nglycan(c):
                                continue
                            out.append({
                                "composition": c,
                                "neutral": masses.neutral_mass(c),
                                "formula": masses.formula_str(c),
                                "name": composition_name(c),
                            })
    return out


def composition_name(c: dict) -> str:
    """카운트로부터 이름 재생성 (엑셀 E열 이름은 틀렸으므로 신뢰 금지)."""
    parts = []
    for key, label in (("HexNAc", "HexNAc"), ("Hex", "Hex"), ("dHex", "Fuc"),
                       ("Neu5Ac", "Neu5Ac"), ("Neu5Gc", "Neu5Gc"), ("Xyl", "Xyl")):
        if c.get(key, 0):
            parts.append(f"{label}{c[key]}")
    return " ".join(parts) if parts else "(empty)"


def adduct_ions(max_charge: int = 3, cations=("H", "Na", "K")):
    """
    전하 1..max_charge 에 대한 모든 adduct 다중집합을 생성.
    반환: [(adduct_tuple, z), ...]
    예: ("H",),1 / ("H","Na"),2 / ("H","H","K"),3 ...
    """
    out = []
    for z in range(1, max_charge + 1):
        for combo in combinations_with_replacement(cations, z):
            out.append((combo, z))
    return out


def adduct_label(adduct) -> str:
    """('H','H','Na') -> '2H+Na' 식 라벨."""
    from collections import Counter
    cnt = Counter(adduct)
    order = ["H", "Na", "K"]
    parts = []
    for cat in order:
        if cnt[cat]:
            parts.append(f"{cnt[cat] if cnt[cat] > 1 else ''}{cat}")
    return "+".join(parts)
