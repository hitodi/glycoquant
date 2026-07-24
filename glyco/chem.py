"""
화학/질량 계산 (설정 구동)
--------------------------
Config(라벨·단당·cation)로부터 글리칸 중성질량, adduct m/z, 진단이온 m/z를
계산하는 Chemistry 객체. masses.py 는 이 객체의 기본(ProA) 인스턴스를 감싼 shim.
"""

# 표준 monoisotopic 원소질량 (중성 원자)
ELEMENT = {
    "C": 12.0, "H": 1.0078250319, "N": 14.0030740052, "O": 15.9949146221,
    "S": 31.97207069, "P": 30.97376151,
    "Na": 22.98976928, "K": 38.96370649,
}
ELECTRON = 0.00054857991


def formula_mass(formula: dict) -> float:
    return sum(ELEMENT[el] * n for el, n in formula.items())


def ppm_error(obs, theo):
    """관측 m/z 의 이론값 대비 ppm 오차 (부호 유지). obs 가 None 이면 None."""
    return (obs - theo) / theo * 1e6 if obs is not None else None


_H2O = formula_mass({"H": 2, "O": 1})


class Chemistry:
    """설정에서 만든 질량 계산기."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.monos = cfg.monosaccharides            # name -> formula
        self.mono_names = tuple(cfg.monosaccharides.keys())
        self.label = cfg.label
        # 잔기질량(= 자유단당 - H2O)
        self.residue = {n: formula_mass(f) - _H2O for n, f in self.monos.items()}
        # 라벨 부착 질량증가 = (자유라벨 - 부착손실) , 보통 ProA: free - O
        self.label_free = formula_mass(self.label.formula)
        self.label_loss = formula_mass(self.label.attach_loss)
        self.label_delta = self.label_free - self.label_loss
        # cation(bare) 질량
        self.cation = {c: ELEMENT[c] - ELECTRON for c in cfg.cations}

    # ---- 중성질량 ----
    def neutral_mass(self, composition: dict) -> float:
        """글리칸 조성 -> 라벨 부착 중성 monoisotopic 질량."""
        res = sum(self.residue[n] * composition.get(n, 0) for n in self.mono_names)
        label_n = composition.get("ProA", composition.get(self.label.name, self.label.count))
        return res + _H2O + self.label_delta * label_n

    def element_formula(self, composition: dict) -> dict:
        tot = {}
        n_sug = 0
        for n in self.mono_names:
            k = composition.get(n, 0)
            n_sug += k
            for el, c in self.monos[n].items():
                tot[el] = tot.get(el, 0) + c * k
        label_n = composition.get("ProA", composition.get(self.label.name, self.label.count))
        for el, c in self.label.formula.items():
            tot[el] = tot.get(el, 0) + c * label_n
        bonds = max(n_sug - 1, 0)
        tot["H"] = tot.get("H", 0) - 2 * bonds
        tot["O"] = tot.get("O", 0) - bonds
        for el, c in self.label.attach_loss.items():
            tot[el] = tot.get(el, 0) - c * label_n
        return {k: v for k, v in tot.items() if v}

    def formula_str(self, composition: dict) -> str:
        f = self.element_formula(composition)
        return "".join(f"{el}{f[el]}" for el in ("C", "H", "N", "O", "S", "P") if f.get(el))

    # ---- 이온 ----
    def ion_mz(self, neutral: float, adduct, z: int) -> float:
        if len(adduct) != z:
            raise ValueError("adduct 개수와 전하수 z 불일치")
        return (neutral + sum(self.cation[a] for a in adduct)) / z

    # ---- 진단이온 ----
    def diagnostic_mz(self, di) -> float:
        """DiagnosticIon -> 단일전하 m/z."""
        if di.kind == "oxonium":      # 잔기합 + proton (B형, 탈수)
            res = sum(self.residue[n] * k for n, k in di.monos.items())
            return res + self.cation["H"]
        elif di.kind == "reducing":   # 라벨 포함 환원말단 단편 + proton (Y형)
            comp = dict(di.monos)
            return self.ion_mz(self.neutral_mass(comp), ("H",), 1)
        raise ValueError(f"알 수 없는 진단이온 kind: {di.kind}")

    def diagnostic_table(self):
        """[(name, m/z), ...]"""
        return [(di.name, self.diagnostic_mz(di)) for di in self.cfg.diagnostic_ions]
