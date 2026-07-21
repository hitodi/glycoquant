"""
진단규칙(diagnostic rule) 파일 로딩 — v2
----------------------------------------
글리칸 공통 단일 검출규칙 + feature 주석. 사용자가 직접 계산한 이론 m/z.
손입력이므로 엄격 검증 + 파싱 결과 echo.

YAML:
  ppm: 5
  ions:                      # 표시·판정에 쓰는 named 코어 이온
    ProA-HexNAc: 441.2708
    HexNAc: 204.0867
    HexNAc+Hex: 366.1395
  accept:                    # 모든 그룹 만족해야 채택(그룹 = named 이온 중 min개 이상)
    - {any: [ProA-HexNAc], min: 1}         # 441 필수
    - {any: [HexNAc, HexNAc+Hex], min: 1}  # 204/366 중 1
  features:                  # 검출 시 주석(채택 무관). mz 리스트면 OR
    - {name: Neu5Ac, mz: [292.1027]}
    - {name: bisecting-GlcNAc, mz: [1009.4823, 1155.5403]}
  group:                     # 구조찾기 그룹핑
    ppm: 5
    split_by_charge: true

⚠️ 441 완화(시알산 손실 우려)는 accept 첫 그룹에 Neu5Ac/Neu5Gc 이온 추가로 처리.
"""

import os
import yaml


class DiagnosticSpec:
    def __init__(self, ions, accept, features, ppm=5.0, group_ppm=5.0, split_by_charge=True):
        self.ions = ions            # {name: m/z}
        self.accept = accept        # [{"any":[name,...], "min":int}]
        self.features = features    # [{"name":str, "mz":[float,...]}]
        self.ppm = ppm
        self.group_ppm = group_ppm
        self.split_by_charge = split_by_charge


def load(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"진단규칙 파일이 없습니다: {path}")
    with open(path, "r", encoding="utf-8") as f:
        d = yaml.safe_load(f) or {}

    ions = {}
    for name, mz in (d.get("ions") or {}).items():
        try:
            ions[str(name)] = float(mz)
        except (TypeError, ValueError):
            raise ValueError(f"이온 '{name}' 값이 숫자가 아닙니다: {mz!r}")

    accept = []
    for grp in (d.get("accept") or []):
        anys = [str(x) for x in (grp.get("any") or [])]
        accept.append({"any": anys, "min": int(grp.get("min", 1))})

    features = []
    for ft in (d.get("features") or []):
        mzs = ft.get("mz")
        mzs = mzs if isinstance(mzs, list) else [mzs]
        try:
            mzs = [float(x) for x in mzs]
        except (TypeError, ValueError):
            raise ValueError(f"feature '{ft.get('name')}' 의 mz 가 숫자가 아닙니다: {ft.get('mz')!r}")
        features.append({"name": str(ft["name"]), "mz": mzs})

    g = d.get("group", {}) or {}
    spec = DiagnosticSpec(
        ions=ions, accept=accept, features=features,
        ppm=float(d.get("ppm", 5.0)),
        group_ppm=float(g.get("ppm", d.get("ppm", 5.0))),
        split_by_charge=bool(g.get("split_by_charge", True)),
    )
    _validate(spec)
    return spec


def _validate(spec):
    if not spec.ions:
        raise ValueError("코어 이온(ions)이 하나도 없습니다.")
    if not spec.accept:
        raise ValueError("채택 규칙(accept)이 없습니다.")
    for grp in spec.accept:
        if not grp["any"]:
            raise ValueError("accept 그룹의 any 가 비었습니다.")
        for n in grp["any"]:
            if n not in spec.ions:
                raise ValueError(f"accept 의 '{n}' 이 ions 에 정의되지 않았습니다.")
        if grp["min"] < 1 or grp["min"] > len(grp["any"]):
            raise ValueError(f"accept min({grp['min']})이 any 개수({len(grp['any'])}) 범위를 벗어남.")


def summary(spec):
    lines = [f"[진단규칙] ppm=±{spec.ppm} | 코어이온 {len(spec.ions)} | feature {len(spec.features)} | 그룹 ppm=±{spec.group_ppm}, charge분리={spec.split_by_charge}"]
    for name, mz in spec.ions.items():
        lines.append(f"   이온 {name}: {mz:.4f}")
    for grp in spec.accept:
        lines.append(f"   채택: [{', '.join(grp['any'])}] 중 {grp['min']}개 이상")
    for ft in spec.features:
        lines.append(f"   feature {ft['name']}: {', '.join(f'{m:.4f}' for m in ft['mz'])}")
    return "\n".join(lines)
