"""
반복 취합 (replicate aggregation)
---------------------------------
여러 파일(반복)의 동정·정량 결과를 글리칸(조성) 단위로 묶어
평균 ± 표준편차 + 검출빈도(n/N)를 계산한다.

- 상대%는 각 파일에서 자기 total 로 정규화된 값이라 파일 간 비교가 가능.
- ⭐결측 처리(재현성 관점): 어떤 반복에서 검출 안 된 글리칸은 그 반복에서 **0%**로 간주.
  → 평균/SD 는 항상 '전체 반복(n_total)' 기준으로 계산. 이래야
    (a) 산발적 검출(1/N)이 SD 로 드러나 낮게 평가되고(false precision 방지),
    (b) 글리칸-레벨 합 == 유형-레벨 값 이 되어 두 시트가 일치한다.
- 검출빈도(n/N)는 별도로 함께 보고(몇 번 실제로 잡혔는지).
- 표준편차는 표본표준편차(ddof=1); n_total<2 면 None(표에 '-').
"""

import numpy as np


def _key(c):
    return (c["HexNAc"], c["Hex"], c["dHex"], c.get("Neu5Ac", 0), c.get("Neu5Gc", 0))


def aggregate(per_file):
    """
    per_file : [(label, results), ...]   results = identify.run() 반환 리스트
    반환: dict(
        files=[label…],
        glycans=[{oxford,name,type,composition,sialylated,fucosylated,
                  pct_by_file{label:pct}, mean,sd,cv,n_detected,n_total,
                  evidence_by_file, evidence}…]  (mean 내림차순),
        types={type: {mean,sd, by_file{label:pct}}},
    )
    """
    files = [label for label, _ in per_file]
    n_total = len(files)

    glymap = {}          # key -> record
    type_by_file = {label: {} for label, _ in per_file}   # label -> {type: pct합}

    for label, results in per_file:
        for r in results:
            k = _key(r["composition"])
            g = glymap.get(k)
            if g is None:
                g = glymap[k] = {
                    "oxford": r.get("oxford", ""),
                    "name": r.get("name", ""),
                    "type": r.get("type", ""),
                    "composition": r["composition"],
                    "sialylated": r.get("sialylated", False),
                    "fucosylated": r.get("fucosylated", False),
                    "pct_by_file": {},
                    "evidence_by_file": {},
                }
            g["pct_by_file"][label] = r.get("relative_pct", 0.0)
            g["evidence_by_file"][label] = r.get("evidence", "MS2")
            type_by_file[label][r["type"]] = type_by_file[label].get(r["type"], 0.0) + r.get("relative_pct", 0.0)

    def _mean_sd(values):
        arr = np.array(values, dtype=float)   # 항상 n_total 길이(결측=0)
        mean = float(arr.mean())
        sd = float(arr.std(ddof=1)) if n_total >= 2 else None
        return mean, sd

    glycans = []
    for k, g in glymap.items():
        # 결측 반복은 0 으로 채운 '전체 반복' 벡터
        vals_full = [g["pct_by_file"].get(f, 0.0) for f in files]
        mean, sd = _mean_sd(vals_full)
        cv = (sd / mean * 100) if (sd is not None and mean) else 0.0
        evs = set(g["evidence_by_file"].values())
        evidence = "MS2" if evs == {"MS2"} else ("+".join(sorted(evs)) if evs else "")
        g.update(mean=mean, sd=sd, cv=cv,
                 n_detected=len(g["pct_by_file"]), n_total=n_total, evidence=evidence)
        glycans.append(g)
    glycans.sort(key=lambda x: x["mean"], reverse=True)

    # 유형별 취합 (글리칸과 동일 정책: 결측=0, 전체 n 기준)
    types = {}
    for t in ("High-mannose", "Hybrid", "Complex"):
        vals = [type_by_file[label].get(t, 0.0) for label in files]
        mean, sd = _mean_sd(vals)
        types[t] = {"mean": mean, "sd": sd,
                    "by_file": {label: type_by_file[label].get(t, 0.0) for label in files}}

    return {"files": files, "glycans": glycans, "types": types}
