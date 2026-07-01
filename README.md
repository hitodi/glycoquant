# glycan_tool — 라벨 글리칸 자동 분석기 (v0.2)

Thermo `.raw`(또는 `.mzML`)에서 **환원말단 표지(ProA 등) 글리칸을 자동으로
동정·정량**해 엑셀로 내보낸다. Xcalibur(윈도우 전용) 없이 macOS/Windows에서 동작.

방법론은 같은 그룹의 논문
*Moon et al., Enzyme and Microbial Technology 196 (2026) 110830*
("Structural, quantitative... N-glycans ... LC-ESI-HCD-MS/MS")을 따른다.

기본 설정은 **ProA N-글리칸**이지만, **`configs/*.yaml` 한 파일만 바꾸면**
다른 라벨(2AB·2AA·free)·단당·adduct·탐색범위로 재사용할 수 있다.

---

## 파이프라인

```
.raw ─(ThermoRawFileParser)→ mzML ─(파싱: MS1 + MS2 단편)→
  설정의 단당·라벨로 조성 조합 생성 + 이론 m/z(H/Na/K, 1~3가)
  → MS2 precursor ppm 매칭으로 동정
  → 진단 oxonium 이온(204.09 HexNAc, 441.27 ProA-HexNAc…)으로 확인
  → MS1 정밀질량 게이트(<5 ppm)
  → 각 adduct XIC '면적' 합산 → 상대량(%)          (논문 §2.5)
  → High-mannose/Hybrid/Complex 분류 + Oxford 명명(FA2, M5…)
  → 엑셀(Glycans / Summary / Adducts / Screening)
```

**Screening 시트** = Xcalibur 에서 스캔마다 손으로 옮기던 그 표를 자동 생성.
각 MS2 스캔의 Scan No · RT · 진단이온(204/366/292/308/441/587) 관측 m/z ·
Precursor m/z · Charge 를 전부 덤프(HexNAc 204 보이는 글리칸 스캔만). `--no-screening` 으로 끔.

### 검증 (tests/)
| 항목 | 기준 | 현재 |
|---|---|---|
| 질량엔진 (계산기 44개 [M+H]+) | < 0.001 ppm | 0.0007 ppm ✅ |
| 논문 진단이온 6종 | < 2 ppm | < 1 ppm ✅ |
| Oxford 명명 (논문 Table 1) | 일치 | 8/8 ✅ |
| 정량 순위 상관 (시트 56조성) | Spearman ≥ 0.57 | +0.60 (area) ✅ |

`pytest` 로 전부 재현 가능(대용량 mzML은 있을 때만 정량 테스트 실행).

---

## 빠른 시작

### 비개발자 (더블클릭)
- **macOS**: `run.command` 더블클릭 → `.raw` 경로 입력(드래그 가능)
- **Windows**: `run.bat` 더블클릭 → `.raw` 경로 입력(드래그 가능)

최초 1회 필요한 패키지·변환기를 자동 설치한다.

### 명령줄
```bash
python -m pip install -r requirements.txt        # 최초 1회
# (윈도우) 변환기:  python setup_thermo.py

python glycan_analyze.py 시료.raw                 # 기본(ProA) 분석
python glycan_analyze.py 시료.raw --ms1-ppm 3     # 더 엄격하게
python glycan_analyze.py 시료.raw -c configs/2ab_nglycan.yaml   # 다른 라벨
```

주요 옵션: `-c/--config`, `--precursor-ppm`, `--ms1-ppm`, `--max-charge`,
`--quant {area,apex}`, `--rt-window`, `--no-ms2`, `--no-diagnostic`,
`--min-intensity`, `--keep-mzml`, `--ms1-first`.

#### `--ms1-first` (시알산 회수 모드, opt-in)
기본 동정은 **MS2 조각화된** 글리칸만 잡는다. 그런데 시알산(특히 Neu5Gc) 글리칸은
Top5 ddMS2 가 잘 안 골라서 **놓치는 경우가 많다**(논문도 이 문제로 MS1 EIC 를 씀).
`--ms1-first` 를 켜면 **MS2 없이도 "같은 RT 에 H·Na·K 등 adduct 가 3개 이상 함께
용출"** 하면 글리칸으로 인정한다(시알산이 알칼리 adduct 를 여럿 만드는 성질을 역이용).

- 검증(논문 PPE 정답 대비): NeuGc 를 **1.0% → 15.7%** 로 회수(논문 14.0%).
- 음성대조군(시알산 없는 시료)에 **위양성을 추가로 만들지 않는 값(noise×70)** 으로 기본 설정.
- ⚠️ 그래도 **정밀도↔회수 트레이드오프**가 있는 실험적 모드다. **시알산이 많은 시료에서만** 권장.
- 결과 엑셀 `Evidence` 열에 각 글리칸이 `MS2` / `MS1` / `MS2+MS1` 중 무엇으로 잡혔는지 표시.

세부 조절: `--ms1-min-adducts`(기본 3), `--ms1-first-ppm`(3), `--ms1-noise-factor`(70).

### 파이썬에서 직접
```python
from glyco import pipeline
results, out_path, cfg = pipeline.analyze("시료.raw", config_path=None)
```

---

## 설정으로 일반화 (`configs/*.yaml`)

라벨·단당·adduct·탐색범위·허용오차·정량방식·진단이온을 전부 YAML로 정의한다.
새 실험은 `proa_nglycan.yaml`을 복사해 수정하면 끝.

```yaml
label:                       # 환원말단 표지
  name: ProA
  formula: {C: 13, H: 21, N: 3, O: 1}
  attach_loss: {O: 1}        # 환원아민화 = -O
monosaccharides:             # 자유 단당 원소조성
  HexNAc: {C: 8, H: 15, N: 1, O: 6}
  Hex:    {C: 6, H: 12, O: 6}
  ...
adducts: {cations: [H, Na, K], max_charge: 3}
search_ranges: {HexNAc: [2,8], Hex: [3,12], ...}
tolerances: {precursor_ppm: 10, ms1_ppm: 5, ms2_ppm: 20}
quantify: {method: area, rt_window_min: 0.5, require_ms2: true}
diagnostic_ions: [...]       # oxonium / reducing 단편
```

`configs/2ab_nglycan.yaml` = 라벨만 2AB로 바꾼 템플릿 예시.

---

## 구조
```
glycan_tool/
  glycan_analyze.py      # CLI (얇은 진입점; 더블클릭 래퍼가 호출)
  setup_thermo.py        # OS별 변환기 자동 설치
  run.command / run.bat  # 비개발자용 더블클릭 래퍼
  pyproject.toml         # 패키징/테스트 설정
  configs/               # *.yaml 설정 (라벨·단당·범위…)
  glyco/                 # 라이브러리
    config.py            #   YAML -> 설정객체
    chem.py              #   질량/이온/진단이온 계산 (설정 구동)
    masses.py            #   chem 의 기본(ProA) shim (하위호환)
    compositions.py      #   조성 조합 생성 + adduct
    nomenclature.py      #   Oxford 명명(FA2, M5…)
    classify.py          #   High-mannose/Hybrid/Complex
    raw.py               #   .raw -> mzML
    mzml_parse.py        #   mzML 파싱 (MS1 + MS2 단편)
    quantify.py          #   XIC apex/area
    diagnostic.py        #   oxonium 확인
    identify.py          #   precursor 매칭 + 확인 + 집계
    pipeline.py          #   오케스트레이션 (analyze())
    report.py            #   엑셀 출력
  tests/                 # 회귀 가드레일(질량/정량/명명/설정)
```

---

## 알아둘 점 / 한계 (v0.2)
- **동정**은 결정론적 m/z + 진단이온 확인이라 신뢰 가능.
- **정량(%)**은 자동 EIC 면적이라 사람이 읽은 값과 절대값은 다를 수 있음
  → 순위·패턴으로 해석(시트와 Spearman ≈ 0.6).
- 진단 oxonium 은 '글리칸 vs 비글리칸'을 거르는 것(노이즈 제거). '조성 vs 조성'
  구분은 못 하므로 동질량 오탐은 `--ms1-ppm` 강화로 줄인다.
- 구조 분류·Oxford 명명은 **조성 기반 근사**(이성질체·안테나 위치 구분 불가).
- **절대정량(pmol)은 범위 밖** — 논문처럼 UPLC 형광 + 표준품 검량선이 필요(MS raw만으론 불가).
- 단일 실행파일 배포: 각 OS에서 `pyinstaller --onefile --add-data "glyco/vendor:glyco/vendor" --add-data "configs:configs" glycan_analyze.py` (크로스컴파일 불가).
