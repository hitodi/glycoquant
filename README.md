# 라벨 N-글리칸 자동 분석기

Thermo `.raw`(또는 `.mzML`)에서 **환원말단 표지(ProA 등) N-글리칸을 자동으로
동정·정량**해 엑셀로 내보내는 CLI 도구. **Xcalibur(윈도우 전용) 없이 macOS/Windows에서** 동작한다.

Xcalibur 에서 스캔 값을 손으로 옮기고, 이론질량을 계산기 시트에 넣고, adduct 강도를
하나씩 찾아 합산하던 **수작업 전 과정을 몇 초로 대체**한다.

> 방법론 근거: *Moon et al., "Structural, quantitative, and functional characterization
> of N-glycans in porcine pancreatin extract, lipase, and α-amylase",
> Enzyme and Microbial Technology 196 (2026) 110830* (LC-ESI-HCD-MS/MS, ProA 표지).

기본 설정은 **ProA N-글리칸**이지만 **`configs/*.yaml` 한 파일만 바꾸면**
다른 라벨(2AB·2AA·free)·단당·adduct·탐색범위로 재사용할 수 있다.

---

## 이게 대체하는 "수작업" (검증됨)

원본 엑셀에서 손으로 하던 작업을, 원본 값과 대조해 커버를 확인했다.

| 손으로 하던 것 | 도구 | 대조 검증 |
|---|---|---|
| **Xcalibur 스캔값 전사** (스크리닝 시트) | `Screening` 시트 자동 | 204 oxonium **12/12 정확 일치**, 글리칸 후보 스캔 자동(전체는 `--screening-all`) |
| **이론 monoisotopic 질량 계산** (계산기) | 질량 엔진 | 계산기 44개 [M+H]⁺ **< 0.001 ppm** |
| **adduct m/z 나열** (H/Na/K × 1~3가) | 자동 생성 | **19종 = 19종 일치** |
| **adduct 강도 찾아 합산** | XIC 면적 자동 | 아래 실데이터 검증 참조 |
| **상대 % / PPM 매칭 / 유형 분류** | 자동 | 유형분포 수작업과 일치 |

---

## 실데이터 검증

도구가 낸 결과를 **사람이 만든 정답(수작업 시트 / 논문 Table)** 과 조성 단위로 대조한 값이다.

**① 알츠하이머 두정엽 (수작업 시트 대조):**
- 공통 글리칸 40개, **Pearson r = 0.94**, 중앙오차 **0.27 %p**, Spearman 0.75
- 상위 글리칸 근접: FA3 21.6 %(수작업 19.8), M5 20.3 %(18.7), FM3 3.61 %(3.53)

**② Porcine pancreatin = 논문 PPE (Table 1 정답 대조):**
- recall 24/32, FA2 **17.4 %**(논문 17.7), 유형분포 HM 30 / Hyb 3 / Cplx 68 %(논문 31/4/64)
- `--ms1-first` 시 Neu5Gc **1.0 % → 15.7 %** 회수(논문 14.0)

**③ 단위 테스트(`pytest`, tests/):**

| 항목 | 기준 | 현재 |
|---|---|---|
| 질량엔진 (계산기 44개 [M+H]⁺) | < 0.001 ppm | 0.0007 ppm ✅ |
| 논문 진단이온 6종 | < 2 ppm | < 1 ppm ✅ |
| Oxford 명명 (논문 Table 1) | 일치 | 8/8 ✅ |
| 정량 순위 상관 (시트 조성, 출하 경로) | Spearman ≥ 0.70 | +0.76 ✅ |

---

## 파이프라인

```
.raw ─(ThermoRawFileParser)→ mzML ─(파싱: MS1 + MS2 단편)→
  설정의 단당·라벨로 조성 조합 생성 + 이론 m/z(H/Na/K, 1~3가)
  → MS2 precursor ppm 매칭으로 동정
  → 진단 oxonium 이온(204.09 HexNAc, 441.27 ProA-HexNAc…)으로 확인
  → MS1 정밀질량 게이트(< 5 ppm)
  → 각 adduct XIC '면적' 합산 → 상대량(%)          (논문 §2.5)
  → High-mannose / Hybrid / Complex 분류 + Oxford 명명(FA2, M5…)
  → 엑셀 4시트: Glycans / Summary / Adducts / Screening
```

---

## 설치

**의존성** (Python ≥ 3.9):
```bash
python -m pip install -r requirements.txt
# numpy, openpyxl, pyteomics, psims, lxml, pyyaml
```

**변환기(ThermoRawFileParser)**: macOS-arm64 빌드는 저장소에 포함(`glyco/vendor/`,
`.gitignore` 됨 — 최초 1회 `python setup_thermo.py` 로 각 OS 빌드 자동 설치).
macOS/Linux 는 `.NET`(dotnet) 런타임 필요, Windows 빌드는 self-contained.

---

## 사용법

### 비개발자 (더블클릭)
- **macOS**: `run.command` 더블클릭 → `.raw` 경로 입력(드래그 가능)
- **Windows**: `run.bat` 더블클릭 → `.raw` 경로 입력(드래그 가능)

최초 1회 필요한 패키지·변환기를 자동 설치한다.

### 명령줄
```bash
python glycan_analyze.py 시료.raw                       # 기본 분석(4시트 + Screening 포함)
python glycan_analyze.py 시료.raw --screening-only       # Xcalibur 스크리닝 전사만 빠르게
python glycan_analyze.py 시료.raw --screening-only --screening-all   # 전체 MS2 스캔 통째로
python glycan_analyze.py 시료.raw -o 결과.xlsx --ms1-ppm 3
python glycan_analyze.py 시료.raw --ms1-first            # 시알산 회수
python glycan_analyze.py 시료.raw -c configs/2ab_nglycan.yaml
```
202
### 파이썬에서
```python
from glyco import pipeline
results, out_path, cfg = pipeline.analyze("시료.raw")
```

### 옵션 전체

| 옵션 | 뜻 | 기본 |
|---|---|---|
| `input` | 입력 `.raw` 또는 `.mzML` (필수) | — |
| `-o, --output` | 출력 `.xlsx` | `입력명_glycans.xlsx` |
| `-c, --config` | 설정 YAML | `configs/proa_nglycan.yaml` |
| `--ms1-ppm` | MS1 정밀질량 게이트(ppm) | 5 |
| `--precursor-ppm` | MS2 precursor 매칭(ppm) | 10 |
| `--max-charge` | 최대 전하수 | 3 |
| `--quant {area,apex}` | 정량 방식(면적/피크높이) | area |
| `--rt-window` | 면적적분 RT 폭(분) | 0.5 |
| `--rt-consistency` | adduct 합산 RT 일치 허용(분), `0`=끔 | 1.0 |
| `--min-intensity` | 이 값 미만 adduct 무시 | 0 |
| `--no-diagnostic` | 진단 oxonium 확인 끔 | (켜짐) |
| `--no-ms2` | MS2 근거 없이도 정량 | (필수) |
| `--ms1-first` | MS1-first 모드(시알산 회수, opt-in) | 꺼짐 |
| `--ms1-min-adducts` / `--ms1-first-ppm` / `--ms1-noise-factor` | MS1-first 세부 | 3 / 3 / 70 |
| `--no-screening` | Screening 시트 생략 | (생성) |
| `--screening-all` | 스크리닝을 필터 없이 **전체 MS2 스캔** 덤프(204 없는 것도 전부) | (204+ 만) |
| `--screening-only` | 동정·정량 없이 Screening 엑셀만 생성 | 꺼짐 |
| `--screening-anchor` / `--screening-ppm` | Screening 포함 기준 진단이온 / 허용오차 | HexNAc,ProA-HexNAc / ms2_ppm |
| `--keep-mzml` | 변환 mzML 보존 | (삭제) |

**상황별 3개만 기억**: 시알산 많은 시료 → `--ms1-first` / 오탐 많음 → `--ms1-ppm 3` / 다른 라벨 → `-c 설정.yaml`.

#### `--ms1-first` (시알산 회수 모드)
기본 동정은 **MS2 조각화된** 글리칸만 잡는다. 시알산(특히 Neu5Gc)은 Top5 ddMS2 가 잘
안 골라 **놓치기 쉽다**(논문도 이 문제로 MS1 EIC 사용). `--ms1-first` 는 MS2 없이도
**"같은 RT 에 H·Na·K 등 adduct 3개 이상 함께 용출"** 하면 글리칸으로 인정한다(시알산이
알칼리 adduct 여럿 만드는 성질을 역이용). 강도 하한은 노이즈 median × 배수(장비 무관 일반화).
- 음성대조군(시알산 없는 시료)에 위양성을 **추가로 안 만드는 값(noise×70)** 으로 기본 설정.
- ⚠️ 정밀도↔회수 트레이드오프가 있는 실험적 모드. 결과 `Evidence` 열에 `MS2`/`MS1`/`MS2+MS1` 표시.

#### 스크리닝 (Xcalibur 수작업 전사 대체)
Xcalibur 에서 스캔마다 Scan No·RT·oxonium m/z·precursor 를 손으로 옮기던 표를 자동 생성한다.
- 기본 분석을 돌리면 결과 엑셀에 **`Screening` 시트가 자동 포함**된다.
- **`--screening-only`**: 동정·정량은 건너뛰고 **스크리닝 엑셀만** 빠르게 만든다(`.raw` 는 MS2 만 변환해 더 빠름). 출력은 `입력명_screening.xlsx`.
- **필터**: 기본은 anchor 진단이온(`HexNAc,ProA-HexNAc`)이 **하나라도 보이는** 글리칸 후보 스캔만.
  - **`--screening-all`**: 필터 없이 **전체 MS2 스캔**을 그대로 덤프(예: 알츠하이머 2080 → 8721 스캔).
  - `--screening-anchor "HexNAc"` / `--screening-ppm 20`: anchor 이온·허용오차 조절(빈 문자열이면 전체).
- **컬럼**: Scan No · RT · 각 진단이온(204/366/292/308/441/587) 관측 m/z · Precursor m/z ·
  Monoisotope/selected m/z · Isolation target · Charge · 비고.

---

## 결과 엑셀
기본 분석은 4시트 — 스크리닝 전용(`--screening-only`) 은 Screening 1시트만.
- **Glycans**: 동정·정량된 글리칸 목록(Oxford명·조성·유형·m/z·RT·Evidence·상대%)
- **Summary**: 유형별(High-man/Hybrid/Complex) + 시알/푸코 집계
- **Adducts**: 글리칸별 검출 adduct 상세
- **Screening**: MS2 스캔별 진단이온 관측 m/z + precursor/monoisotope/isolation/charge (수작업 전사 대체)

---

## 설정으로 일반화 (`configs/*.yaml`)

라벨·단당·adduct·탐색범위·허용오차·정량방식·진단이온을 전부 YAML로 정의한다.
새 실험은 `proa_nglycan.yaml`을 복사해 수정하면 끝. `configs/2ab_nglycan.yaml` = 2AB 라벨 예시.

```yaml
label:                       # 환원말단 표지
  name: ProA
  formula: {C: 13, H: 21, N: 3, O: 1}
  attach_loss: {O: 1}        # 환원아민화 = -O
monosaccharides:
  HexNAc: {C: 8, H: 15, N: 1, O: 6}
  Hex:    {C: 6, H: 12, O: 6}
  # ...
adducts: {cations: [H, Na, K], max_charge: 3}
search_ranges: {HexNAc: [2, 8], Hex: [3, 12]}   # ...
tolerances: {precursor_ppm: 10, ms1_ppm: 5, ms2_ppm: 20}
quantify: {method: area, rt_window_min: 0.5, require_ms2: true}
diagnostic_ions: [...]       # oxonium / reducing 단편
```

---

## 구조
```
glycan_analyze.py      # CLI 진입점 (더블클릭 래퍼가 호출)
setup_thermo.py        # OS별 변환기 자동 설치
run.command / run.bat  # 비개발자용 더블클릭 래퍼
pyproject.toml         # 패키징/테스트 설정
configs/               # *.yaml 설정 (proa_nglycan, 2ab_nglycan)
glyco/                 # 라이브러리
  config.py            #   YAML → 설정객체
  chem.py              #   질량/이온/진단이온 (설정 구동)
  masses.py            #   chem 의 기본(ProA) shim
  compositions.py      #   조성 조합 생성 + adduct + 타당성
  nomenclature.py      #   Oxford 명명(FA2, M5…)
  classify.py          #   유형 분류
  raw.py               #   .raw → mzML
  mzml_parse.py        #   mzML 파싱 (MS1 + MS2 단편)
  quantify.py          #   XIC apex/area
  diagnostic.py        #   oxonium 확인
  identify.py          #   precursor 매칭 + 확인 + 집계
  pipeline.py          #   오케스트레이션 analyze()
  report.py            #   엑셀 출력
scripts/               # 검증 하네스(eval_ppe: 논문정답 / fp_probe: 위양성)
tests/                 # 회귀 가드레일(질량/정량/명명/설정)
```

---

## 알아둘 점 / 한계
- **동정**은 결정론적 m/z + 진단이온 확인이라 신뢰 가능.
- **정량(%)**은 자동 EIC 면적이라 사람이 읽은 값과 절대값은 다를 수 있음
  → 순위·패턴으로 해석(수작업과 Pearson ≈ 0.94, 중앙오차 ~0.3 %p).
- 진단 oxonium 은 '글리칸 vs 비글리칸' 필터(노이즈 제거). 동질량 오탐은 `--ms1-ppm` 강화로 감소.
- 구조 분류·Oxford 명명은 **조성 기반 근사**(이성질체·안테나 위치 구분 불가).
- **절대정량(pmol)은 범위 밖** — 논문처럼 UPLC 형광 + 표준품 검량선 필요(MS raw만으론 불가).
- 단일 실행파일 배포: 각 OS에서
  `pyinstaller --onefile --add-data "glyco/vendor:glyco/vendor" --add-data "configs:configs" glycan_analyze.py`
  (PyInstaller 는 크로스컴파일 불가 → 대상 OS 에서 빌드).

---

## 인용 / 라이선스
- 방법론 인용: Moon et al., *Enzyme and Microbial Technology* 196 (2026) 110830.
- 번들: [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) (자체 라이선스).
- 라이선스: [MIT](LICENSE)
