# 사용방법

이 문서는 설치, 실행 예시, 옵션, `--targets` 진단 스크리닝 YAML 작성법을 정리한다.

---

## 설치

Python 3.9 이상이 필요하다.

```bash
python -m pip install -r requirements.txt
```

Thermo `.raw` 입력은 ThermoRawFileParser로 `.mzML`로 변환된다. 최초 1회 다음 명령으로
OS별 변환기를 준비한다.

```bash
python setup_thermo.py
```

macOS/Linux에서는 `.NET` 런타임이 필요하다. Windows 빌드는 self-contained 형태다.

---

## 빠른 실행

```bash
python glycan_analyze.py 시료.raw
python glycan_analyze.py 시료.mzML
python glycan_analyze.py 반복폴더/
```

기본 출력:

- 파일 입력: `시료_glycans.xlsx`
- 폴더 입력: `반복폴더/glycoquant_results/` 아래 파일별 결과와 `aggregated.xlsx`

출력 경로를 직접 지정하려면:

```bash
python glycan_analyze.py 시료.raw -o 결과.xlsx
python glycan_analyze.py 반복폴더/ -o 결과폴더/
```

---

## 비개발자 실행

- macOS: `run.command` 더블클릭 후 `.raw` 경로 입력 또는 파일 드래그
- Windows: `run.bat` 더블클릭 후 `.raw` 경로 입력 또는 파일 드래그

최초 실행 시 필요한 패키지와 변환기를 자동 설치한다.

---

## 자주 쓰는 명령

```bash
# 파일 1개 분석
python glycan_analyze.py 시료.raw

# 폴더의 .raw/.mzML 여러 개를 반복 실험으로 취합
python glycan_analyze.py 반복폴더/ -o 결과폴더/

# Xcalibur 수작업 전사용 Screening 엑셀만 빠르게 생성
python glycan_analyze.py 시료.raw --screening-only

# 전체 MS2 스캔을 필터 없이 Screening에 덤프
python glycan_analyze.py 시료.raw --screening-only --screening-all

# MS1 정밀질량 게이트를 3 ppm으로 강화
python glycan_analyze.py 시료.raw --ms1-ppm 3

# 시알산 글리칸 회수용 MS1-first 모드
python glycan_analyze.py 시료.raw --ms1-first

# 2AB 설정으로 분석
python glycan_analyze.py 시료.raw -c configs/2ab_nglycan.yaml

# 진단이온 규칙 기반 스크리닝
python glycan_analyze.py 시료.raw --targets configs/diagnostics_proa.yaml
```

상황별로 3개만 기억하면 된다.

- 시알산 많은 시료: `--ms1-first`
- 오탐이 많음: `--ms1-ppm 3`
- 다른 라벨: `-c 설정.yaml`

---

## 진단 스크리닝: `--targets`

`--targets`는 연구실에서 직접 계산한 진단이온 이론 m/z와 채택 규칙을 YAML로 넣어
MS2 스캔을 선별하는 모드다. 일반 동정·정량을 하지 않고, 스크리닝용 엑셀을 만든다.

```bash
python glycan_analyze.py 시료.raw --targets configs/diagnostics_proa.yaml
```

출력:

```text
시료_screening.xlsx
```

출력 경로 지정:

```bash
python glycan_analyze.py 시료.raw \
  --targets configs/diagnostics_proa.yaml \
  -o 결과.xlsx
```

실행 시 ppm만 덮어쓰기:

```bash
python glycan_analyze.py 시료.raw \
  --targets configs/diagnostics_proa.yaml \
  --targets-ppm 10 \
  --group-ppm 5
```

출력 엑셀은 2개 시트로 구성된다.

- `Screening`: 채택된 MS2 스캔, RT, 코어이온 관측 m/z, ppm 오차, precursor, monoisotope, charge, Features
- `구조찾기`: 채택 스캔을 monoisotope ±ppm, charge별로 그룹핑하고 `SC no.`로 Screening 행을 역참조

monoisotope 열과 그룹키는 Thermo `Monoisotopic M/Z`가 있으면 그 값을 쓰고, 없으면 precursor를 쓴다.

---

## `diagnostics_proa.yaml` 행별 설명

대상 파일: `configs/diagnostics_proa.yaml`

```yaml
ppm: 5
```

진단이온 m/z를 찾을 때 허용할 오차다. 예를 들어 이론값 `441.2708`을 ±5 ppm 안에서 찾는다.

```yaml
ions:
  ProA-HexNAc: 441.2708
  HexNAc: 204.0867
  HexNAc+Hex: 366.1395
```

`ions`는 **채택 판정에 직접 쓰는 핵심 진단이온 목록**이다.

- `ProA-HexNAc`: ProA 라벨이 붙은 HexNAc 조각, 441 이온
- `HexNAc`: 일반적인 N-glycan oxonium ion, 204 이온
- `HexNAc+Hex`: HexNAc-Hex 조각, 366 이온

여기에 정의된 이름만 아래 `accept`에서 사용할 수 있다.

```yaml
accept:
  - {any: [ProA-HexNAc], min: 1}
  - {any: [HexNAc, HexNAc+Hex], min: 1}
```

`accept`는 **어떤 MS2 스캔을 채택할지 정하는 규칙**이다. 여러 줄이 있으면 모든 줄을 만족해야 채택된다.

```yaml
- {any: [ProA-HexNAc], min: 1}
```

`ProA-HexNAc`가 반드시 1개 이상 보여야 한다는 뜻이다. 즉 441 이온 필수다.

```yaml
- {any: [HexNAc, HexNAc+Hex], min: 1}
```

`HexNAc` 또는 `HexNAc+Hex` 중 하나 이상이 보여야 한다는 뜻이다. 즉 204 또는 366 중 하나가 있으면 통과한다.

현재 기본 규칙은 다음과 같다.

```text
441 필수 AND (204 OR 366)
```

```yaml
features:
  - {name: "Neu5Ac(시알산·인간)",   mz: [292.1027]}
  - {name: "Neu5Gc(시알산·비인간)", mz: [308.0976]}
  - {name: "bisecting-GlcNAc",      mz: [1009.4823, 1155.5403]}
  - {name: "core-fucose",           mz: [587.3281]}
  - {name: "terminal-fucose",       mz: [512.1974]}
  - {name: "phospho(Man)",          mz: [243.0259]}
  - {name: "acetyl(Neu5Ac)",        mz: [334.1127]}
  - {name: "acetyl(Neu5Gc)",        mz: [350.1076]}
  - {name: "sulfation(HexNAc)",     mz: [284.0430]}
```

`features`는 **채택 여부에는 영향을 주지 않는 주석용 이온**이다. 해당 m/z가 보이면 결과 엑셀의
`Features` 열에 이름이 붙는다.

예를 들어 어떤 스캔이 441 + 204 조건을 만족해서 채택됐고, 동시에 `292.1027`도 보이면
`Features` 열에 `Neu5Ac(시알산·인간)`이 표시된다.

`mz`에 값이 여러 개 있으면 OR 조건이다.

```yaml
- {name: "bisecting-GlcNAc", mz: [1009.4823, 1155.5403]}
```

이 경우 `1009.4823` 또는 `1155.5403` 중 하나라도 보이면 `bisecting-GlcNAc`로 주석된다.

```yaml
group:
  ppm: 5
  split_by_charge: true
```

`group`은 엑셀의 `구조찾기` 시트에서 스캔들을 묶는 기준이다.

```yaml
ppm: 5
```

monoisotope m/z가 ±5 ppm 안에 있으면 같은 그룹으로 묶는다.

```yaml
split_by_charge: true
```

전하가 다르면 같은 m/z 근처라도 다른 그룹으로 나눈다. 즉 `2+`와 `3+`는 별도 그룹이다.

전체 동작은 다음과 같다.

```text
1. 모든 MS2 스캔을 읽는다.
2. 각 스캔에서 441, 204, 366 이온을 ±5 ppm으로 찾는다.
3. 441이 있고, 204 또는 366이 있으면 채택한다.
4. 채택된 스캔에서 시알산/푸코스/bisecting 등 feature 주석용 이온을 추가로 찾는다.
5. 채택 스캔을 monoisotope ±5 ppm, charge별로 묶어 구조찾기 시트에 정리한다.
```

---

## 시알산 손실 완화 규칙

기본 규칙은 441을 필수로 보기 때문에 일부 시알산 글리칸이 빠질 수 있다. 완화하려면
`ions`에 `Neu5Ac`, `Neu5Gc`를 추가하고 첫 번째 `accept` 그룹에 함께 넣는다.

```yaml
ions:
  ProA-HexNAc: 441.2708
  HexNAc: 204.0867
  HexNAc+Hex: 366.1395
  Neu5Ac: 292.1027
  Neu5Gc: 308.0976

accept:
  - {any: [ProA-HexNAc, Neu5Ac, Neu5Gc], min: 1}
  - {any: [HexNAc, HexNAc+Hex], min: 1}
```

이 규칙은 다음 의미다.

```text
(441 OR 292 OR 308) AND (204 OR 366)
```

회수율은 올라갈 수 있지만 오탐도 늘 수 있으므로, 기본 규칙 결과와 비교해서 쓰는 것이 좋다.

---

## 폴더 배치와 반복 취합

입력에 파일 대신 폴더를 주면, 그 안의 `.raw`/`.mzML` 전부를 각각 분석하고 반복 실험으로 취합한다.

```bash
python glycan_analyze.py 반복폴더/ -o 결과폴더/
```

출력:

- 개별 결과: `결과폴더/<파일명>_glycans.xlsx`
- 취합 결과: `결과폴더/aggregated.xlsx`

`aggregated.xlsx` 시트:

- `Aggregated`: 글리칸별 평균 %, SD, CV %, 검출빈도 n/N, 반복별 %
- `Type summary`: 유형별 평균 ± SD

결측 처리는 재현성 관점으로 한다. 어떤 반복에서 안 잡힌 글리칸은 그 반복을 0%로 간주한다.
산발 검출(예: 1/3)은 SD가 커져 자동으로 신뢰도가 낮게 드러난다.

현재는 **한 폴더 = 한 반복그룹**이다. 여러 시료군을 비교하려면 시료군별로 폴더를 나눠 실행한다.

---

## `--ms1-first`

기본 동정은 MS2 조각화된 글리칸만 잡는다. 시알산, 특히 Neu5Gc는 TopN ddMS2에서 선택되지 않아
놓치기 쉽다. `--ms1-first`는 MS2 없이도 같은 RT에 H/Na/K 등 adduct 3개 이상이 함께 용출되면
글리칸으로 인정한다.

```bash
python glycan_analyze.py 시료.raw --ms1-first
```

세부 옵션:

```bash
python glycan_analyze.py 시료.raw \
  --ms1-first \
  --ms1-min-adducts 3 \
  --ms1-first-ppm 3 \
  --ms1-noise-factor 70
```

정밀도와 회수율 사이의 트레이드오프가 있으므로 시알산이 많은 시료에서 우선 사용한다.
결과 `Evidence` 열에 `MS2`, `MS1`, `MS2+MS1`로 근거가 표시된다.

---

## 기본 Screening 시트

일반 분석을 돌리면 결과 엑셀에 `Screening` 시트가 자동 포함된다.

```bash
python glycan_analyze.py 시료.raw
```

스크리닝 엑셀만 만들려면:

```bash
python glycan_analyze.py 시료.raw --screening-only
```

기본 필터는 anchor 진단이온(`HexNAc,ProA-HexNAc`) 중 하나라도 보이는 MS2 스캔이다.
전체 MS2 스캔을 모두 내보내려면:

```bash
python glycan_analyze.py 시료.raw --screening-only --screening-all
```

anchor 이온과 ppm을 직접 지정하려면:

```bash
python glycan_analyze.py 시료.raw \
  --screening-only \
  --screening-anchor "HexNAc" \
  --screening-ppm 20
```

---

## 분석 설정 YAML

일반 동정·정량 설정은 `configs/*.yaml`로 바꾼다. 새 라벨이나 탐색 범위가 필요하면
`configs/proa_nglycan.yaml`을 복사해서 수정한다. `configs/2ab_nglycan.yaml`은 2AB 라벨 예시다.

```yaml
label:
  name: ProA
  formula: {C: 13, H: 21, N: 3, O: 1}
  attach_loss: {O: 1}

monosaccharides:
  HexNAc: {C: 8, H: 15, N: 1, O: 6}
  Hex:    {C: 6, H: 12, O: 6}

adducts:
  cations: [H, Na, K]
  max_charge: 3

search_ranges:
  HexNAc: [2, 8]
  Hex: [3, 12]

tolerances:
  precursor_ppm: 10
  ms1_ppm: 5
  ms2_ppm: 20

quantify:
  method: area
  rt_window_min: 0.5
  require_ms2: true

diagnostic_ions:
  # oxonium / reducing-end fragment definitions
```

자주 바꾸는 항목:

- `label`: 환원말단 표지 종류와 화학식
- `monosaccharides`: 사용할 단당 조성
- `adducts`: 탐색할 adduct와 최대 전하수
- `search_ranges`: 후보 조성 생성 범위
- `tolerances`: precursor, MS1, MS2 ppm 허용오차
- `quantify`: 면적/apex 정량 방식과 RT 윈도우
- `diagnostic_ions`: 일반 분석에서 확인할 oxonium/reducing-end 단편

---

## 옵션 전체

| 옵션 | 뜻 | 기본 |
|---|---|---|
| `input` | 입력 `.raw`/`.mzML` 또는 폴더 | 필수 |
| `-o, --output` | 출력 `.xlsx` 또는 폴더 배치 출력 폴더 | `입력명_glycans.xlsx` |
| `-c, --config` | 분석 설정 YAML | `configs/proa_nglycan.yaml` |
| `--ms1-ppm` | MS1 정밀질량 게이트(ppm) | 5 |
| `--precursor-ppm` | MS2 precursor 매칭(ppm) | 10 |
| `--max-charge` | 최대 전하수 | 3 |
| `--quant {area,apex}` | 정량 방식 | `area` |
| `--rt-window` | 면적적분 RT 폭(분) | 0.5 |
| `--rt-consistency` | adduct 합산 RT 일치 허용(분), `0`은 끔 | 1.0 |
| `--min-intensity` | 이 값 미만 adduct 무시 | 0 |
| `--no-diagnostic` | 진단 oxonium 확인 끄기 | 켜짐 |
| `--no-ms2` | MS2 근거 없이도 정량 | 꺼짐 |
| `--ms1-first` | MS1-first 모드 | 꺼짐 |
| `--ms1-min-adducts` | MS1-first 공존 adduct 최소 수 | 3 |
| `--ms1-first-ppm` | MS1-first의 MS1 허용오차(ppm) | 3 |
| `--ms1-noise-factor` | MS1-first 강도하한 = 노이즈 median × 값 | 70 |
| `--no-screening` | 일반 분석 결과에서 Screening 시트 생략 | 생성 |
| `--screening-all` | 기본 Screening을 전체 MS2 스캔으로 덤프 | 꺼짐 |
| `--screening-only` | 동정·정량 없이 Screening 엑셀만 생성 | 꺼짐 |
| `--screening-anchor` | Screening 포함 기준 진단이온 | `HexNAc,ProA-HexNAc` |
| `--screening-ppm` | Screening 진단이온 허용오차 | 설정 `ms2_ppm` |
| `--targets` | 진단규칙 YAML 기반 스크리닝 모드 | 없음 |
| `--targets-ppm` | `--targets` 이온 매칭 ppm 덮어쓰기 | 파일값 |
| `--group-ppm` | `--targets` 구조찾기 그룹 ppm 덮어쓰기 | 파일값 |
| `--keep-mzml` | 변환 mzML 보존 | 삭제 |

---

## 파이썬에서 사용

```python
from glyco import pipeline

results, out, cfg = pipeline.analyze("시료.raw")
per_file, agg, agg_path = pipeline.batch("반복폴더/")
rows, out, spec = pipeline.targeted_screening(
    "시료.raw",
    "configs/diagnostics_proa.yaml",
)
```
