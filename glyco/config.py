"""
설정 로딩 (YAML -> dataclass)
-----------------------------
configs/*.yaml 을 읽어 검증된 설정 객체로 변환한다.
이 파일 하나만 바꾸면 라벨(ProA->2AB 등)·단당·adduct·범위·허용오차를 교체할 수 있다.
"""

import os
from dataclasses import dataclass, field

import yaml

CONFIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "configs")
DEFAULT_CONFIG = os.path.join(CONFIG_DIR, "proa_nglycan.yaml")


@dataclass
class Label:
    name: str
    formula: dict
    attach_loss: dict = field(default_factory=dict)
    count: int = 1


@dataclass
class DiagnosticIon:
    name: str
    kind: str          # 'oxonium' | 'reducing'
    monos: dict


@dataclass
class Config:
    name: str
    label: Label
    monosaccharides: dict          # name -> element formula
    cations: list                  # ['H','Na','K']
    max_charge: int
    search_ranges: dict            # name -> (min,max)
    plausibility: dict
    tolerances: dict
    quantify: dict
    diagnostic_ions: list          # [DiagnosticIon]
    require_diagnostic: list

    # 편의 접근자
    @property
    def precursor_ppm(self): return self.tolerances.get("precursor_ppm", 10.0)
    @property
    def ms1_ppm(self): return self.tolerances.get("ms1_ppm", 5.0)
    @property
    def ms2_ppm(self): return self.tolerances.get("ms2_ppm", 20.0)


def load(path: str = None) -> Config:
    path = path or DEFAULT_CONFIG
    with open(path, "r", encoding="utf-8") as f:
        d = yaml.safe_load(f)

    lab = d["label"]
    label = Label(name=lab["name"], formula=dict(lab["formula"]),
                  attach_loss=dict(lab.get("attach_loss", {})), count=int(lab.get("count", 1)))

    ranges = {k: (int(v[0]), int(v[1])) for k, v in d["search_ranges"].items()}
    diag = [DiagnosticIon(name=x["name"], kind=x["kind"], monos=dict(x["monos"]))
            for x in d.get("diagnostic_ions", [])]

    cfg = Config(
        name=d.get("name", os.path.basename(path)),
        label=label,
        monosaccharides={k: dict(v) for k, v in d["monosaccharides"].items()},
        cations=list(d["adducts"]["cations"]),
        max_charge=int(d["adducts"].get("max_charge", 3)),
        search_ranges=ranges,
        plausibility=dict(d.get("plausibility", {})),
        tolerances=dict(d.get("tolerances", {})),
        quantify=dict(d.get("quantify", {})),
        diagnostic_ions=diag,
        require_diagnostic=list(d.get("require_diagnostic", [])),
    )
    _validate(cfg)
    return cfg


def _validate(cfg: Config):
    assert cfg.monosaccharides, "monosaccharides 가 비어 있습니다."
    for c in cfg.cations:
        assert c in ("H", "Na", "K"), f"지원하지 않는 cation: {c} (H/Na/K)"
    for name in cfg.search_ranges:
        assert name in cfg.monosaccharides, f"search_ranges 의 '{name}' 가 monosaccharides 에 없음"
    for di in cfg.diagnostic_ions:
        assert di.kind in ("oxonium", "reducing"), f"diagnostic kind 오류: {di.kind}"
