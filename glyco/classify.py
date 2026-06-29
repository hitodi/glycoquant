"""
구조 분류 (조성 기반 휴리스틱)
------------------------------
구조 정보 없이 조성만으로 분류하므로 근사치다. (보고용 시트 분류 체계에 맞춤)
  - High-mannose : HexNAc==2 (코어만), 시알/푸코 없음, Hex 다수
  - Hybrid       : HexNAc==3
  - Complex      : HexNAc>=4
부가 특성: 시알릴화/푸코실화 여부와 개수.
"""


def classify(c: dict) -> dict:
    hexnac = c["HexNAc"]
    sia = c["Neu5Ac"] + c["Neu5Gc"]
    fuc = c["dHex"]

    if hexnac <= 2:
        typ = "High-mannose"
    elif hexnac == 3:
        typ = "Hybrid"
    else:
        typ = "Complex"

    return {
        "type": typ,
        "sialylated": sia > 0,
        "fucosylated": fuc > 0,
        "n_sialic": sia,
        "n_fucose": fuc,
        "n_neu5ac": c["Neu5Ac"],
        "n_neu5gc": c["Neu5Gc"],
    }
