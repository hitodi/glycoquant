# glycan_tool — ProA 표지 N-글리칸 분석기

Thermo `.raw`(또는 `.mzML`) 파일에서 ProA로 표지한 N-글리칸을
**자동으로 동정(identification)하고 정량(quantification)** 해서 엑셀로 내보냅니다.
Xcalibur(윈도우 전용) 없이 macOS/Windows에서 바로 돌아갑니다.

기존 `Alzhemer of brain.xlsx`에서 손으로 하던
"스크리닝 → 계산기 대조 → 상대량 집계" 과정을 대체합니다.

---

## 무엇을 하나

```
.raw  ─(ThermoRawFileParser)→  mzML  ─(파싱)→  MS1 / MS2 precursor
     → 조성 조합 생성 + 이론 m/z(H/Na/K, 1~3가) 계산
     → MS2 precursor 를 ppm 매칭해 글리칸 동정
     → 각 adduct XIC 피크높이 합산 → 상대량(%)
     → High-mannose / Hybrid / Complex 분류
     → 결과 엑셀(Glycans / Summary / Adducts 3시트)
```

질량 모델은 기존 `계산기` 시트와 동일하게 검증됨(44개 [M+H]+ 재현 오차 < 0.001 ppm).

---

## 빠른 시작

### 비개발자 (더블클릭)
- **macOS**: `run.command` 더블클릭 → `.raw` 경로 입력(드래그도 됨)
- **Windows**: `run.bat` 더블클릭 → `.raw` 경로 입력(드래그도 됨)

최초 1회 자동으로 필요한 패키지와 변환기를 설치합니다.

### 명령줄
```bash
# 의존성(최초 1회)
python -m pip install -r requirements.txt
# (윈도우만) 변환기 받기:  python setup_thermo.py

# 기본 실행 — 입력옆에 같은이름_glycans.xlsx 생성
python glycan_analyze.py 시료.raw

# 옵션 예시
python glycan_analyze.py 시료.raw -o 결과.xlsx --ppm 8 --max-charge 3
```

주요 옵션:

| 옵션 | 뜻 | 기본 |
|---|---|---|
| `-o, --output` | 출력 엑셀 경로 | `입력명_glycans.xlsx` |
| `--ppm` | 질량 허용오차(ppm) | 10 |
| `--max-charge` | 최대 전하수 | 3 |
| `--no-ms2` | MS2 근거 없이도 정량 | 꺼짐(=MS2 확인된 것만) |
| `--min-intensity` | 이 강도 미만 adduct 무시 | 0 |
| `--hexnac/--hex/--fuc/--neu5ac/--neu5gc` | 조성 탐색 범위 `MIN MAX` | 내장 기본값 |
| `--xyl` | Xyl 포함 탐색 | 꺼짐 |
| `--keep-mzml` | 변환 mzML 보존 | 꺼짐 |

---

## 폴더 구조
```
glycan_tool/
  glycan_analyze.py      # CLI 진입점
  setup_thermo.py        # OS별 변환기 자동 설치
  run.command / run.bat  # 더블클릭 실행 래퍼
  requirements.txt
  glyco/
    masses.py            # 질량 엔진(검증됨)
    compositions.py      # 조성 조합 생성 + adduct
    raw.py               # .raw -> mzML
    mzml_parse.py        # mzML 파싱
    quantify.py          # XIC 정량(벡터화)
    identify.py          # MS2 매칭 + 결과 조립
    classify.py          # 구조 유형 분류
    report.py            # 엑셀 출력
    vendor/              # ThermoRawFileParser 번들(.gitignore)
```

---

## 알아둘 점 / 한계 (v0.1)

- **동정**은 결정론적 m/z 매칭이라 신뢰 가능. **정량(%)**은 자동 XIC 값이라
  사람이 Xcalibur에서 읽은 NL 값과 **절대값이 다를 수 있음** → 순위·패턴으로 해석할 것.
- MS2 precursor m/z는 **isolation 타깃값**이라 MS1 정밀질량보다 오차가 큼.
  `--ppm 3`처럼 너무 좁히면 진짜 글리칸을 놓침(기본 10 권장).
- 조성을 조합으로 생성하므로, 넓은 범위 + 느슨한 ppm에서는 **오탐(특히 Complex)** 이
  늘 수 있음. `--ppm` 을 줄이거나 `--min-intensity` 로 약한 신호를 거르면 개선됨.
- 구조 분류는 조성 기반 **휴리스틱**(이성질체 구분 불가).
- 윈도우 배포 시 `setup_thermo.py`가 win 변환기를 받습니다(.NET 포함 self-contained).

### 비개발자에게 단일 실행파일로 배포 (선택)
파이썬조차 없는 PC용으로는 각 OS에서 PyInstaller로 빌드:
```bash
python -m pip install pyinstaller
pyinstaller --onefile --add-data "glyco/vendor:glyco/vendor" glycan_analyze.py
# 윈도우는 구분자 ; :  --add-data "glyco/vendor;glyco/vendor"
```
(PyInstaller는 크로스 컴파일이 안 되므로 윈도우 exe는 윈도우에서 빌드해야 함)
