"""
Oxford 명명법 (조성 기반 근사)
------------------------------
논문 §2.11 Oxford notation: A(안테나 GlcNAc), F(Fuc), G(Gal), M(Man),
P(phosphate), S(Neu5Ac), Sg(Neu5Gc).

⚠️ 조성만으로는 구조(이성질체·안테나 위치)를 알 수 없어 '근사'다.
   trimannosyl 코어(Man3GlcNAc2) 가정으로:
     안테나 A = HexNAc - 2,  갈락토스 G = Hex - 3  (코어 너머 Hex)
   논문 Table 1 의 FA2 / FA2G2S2 / M5 / FA2G1Sg1 등과 일치 확인됨.
"""


def oxford(c: dict) -> str:
    hexnac = c.get("HexNAc", 0)
    hexn = c.get("Hex", 0)
    fuc = c.get("dHex", 0)
    ac = c.get("Neu5Ac", 0)
    gc = c.get("Neu5Gc", 0)

    # high-mannose: 코어 HexNAc2 만, Man = Hex
    if hexnac == 2 and ac == 0 and gc == 0:
        return ("F" if fuc else "") + f"M{hexn}"

    antennae = max(hexnac - 2, 0)
    gal = max(hexn - 3, 0)

    name = ""
    if fuc:
        name += "F" if fuc == 1 else f"F{fuc}"
    name += f"A{antennae}"
    if gal:
        name += f"G{gal}"
    if ac:
        name += f"S{ac}"
    if gc:
        name += f"Sg{gc}"
    return name or f"HexNAc{hexnac}Hex{hexn}"
