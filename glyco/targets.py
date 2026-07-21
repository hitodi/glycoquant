"""
타깃(사용자 지정 진단이온) 파일 로딩
------------------------------------
사용자가 직접 계산한 글리칸별 진단이온 이론 m/z 목록을 읽는다.
손으로 작성하는 입력이므로 **엄격히 검증하고 파싱 결과를 되돌려 출력**한다.

YAML 형식:
  ppm: 20            # (선택) 이온 매칭 허용오차
  min_hits: 2        # (선택) 채택 임계 K
  precursor_floor: 0.1  # (선택) precursor 그룹 버림 단위
  glycans:
    - name: GlycanA
      ions: [204.0867, 366.1395, 528.19]
    - name: GlycanB
      ions: [292.1027, 657.2349, 946.34]

CSV 형식(헤더 없음): 각 행 = name,ion1,ion2,...
"""

import os
import yaml


class TargetSpec:
    def __init__(self, glycans, ppm=20.0, min_hits=2, precursor_floor=0.1):
        self.glycans = glycans          # [{"name","ions":[float,...]}]
        self.ppm = ppm
        self.min_hits = min_hits
        self.precursor_floor = precursor_floor


def load(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"타깃 파일이 없습니다: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".yaml", ".yml"):
        spec = _load_yaml(path)
    elif ext == ".csv":
        spec = _load_csv(path)
    else:
        raise ValueError(f"지원 형식은 .yaml/.yml/.csv 입니다(입력: {ext})")
    _validate(spec)
    return spec


def _num_list(name, raw):
    ions = []
    for v in raw:
        try:
            ions.append(float(v))
        except (TypeError, ValueError):
            raise ValueError(f"글리칸 '{name}' 의 이온값이 숫자가 아닙니다: {v!r}")
    return ions


def _load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        d = yaml.safe_load(f) or {}
    if "glycans" not in d:
        raise ValueError("YAML 에 'glycans' 키가 없습니다.")
    glycans = []
    for i, g in enumerate(d["glycans"]):
        if "name" not in g or "ions" not in g:
            raise ValueError(f"{i+1}번째 글리칸에 name/ions 가 없습니다: {g!r}")
        glycans.append({"name": str(g["name"]), "ions": _num_list(g["name"], g["ions"])})
    return TargetSpec(glycans,
                      ppm=float(d.get("ppm", 20.0)),
                      min_hits=int(d.get("min_hits", 2)),
                      precursor_floor=float(d.get("precursor_floor", 0.1)))


def _load_csv(path):
    glycans = []
    with open(path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                raise ValueError(f"{ln}행: 'name,ion1,...' 형식이 아닙니다: {line!r}")
            glycans.append({"name": parts[0], "ions": _num_list(parts[0], parts[1:])})
    return TargetSpec(glycans)


def _validate(spec):
    if not spec.glycans:
        raise ValueError("글리칸이 하나도 없습니다.")
    if spec.min_hits < 1:
        raise ValueError(f"min_hits 는 1 이상이어야 합니다(입력 {spec.min_hits}).")
    for g in spec.glycans:
        if not g["ions"]:
            raise ValueError(f"글리칸 '{g['name']}' 의 이온 목록이 비었습니다.")
        if len(g["ions"]) < spec.min_hits:
            raise ValueError(
                f"글리칸 '{g['name']}' 이온 {len(g['ions'])}개 < min_hits {spec.min_hits} "
                f"(이 글리칸은 절대 채택될 수 없음).")


def summary(spec):
    """파싱 결과를 사람이 눈으로 검증할 수 있게 문자열로."""
    lines = [f"[타깃] 글리칸 {len(spec.glycans)}개 | ppm={spec.ppm} | min_hits(K)={spec.min_hits} | precursor_floor={spec.precursor_floor}"]
    for g in spec.glycans:
        ions = ", ".join(f"{x:.4f}" for x in g["ions"])
        lines.append(f"   - {g['name']}: {len(g['ions'])}개 [{ions}]")
    return "\n".join(lines)
