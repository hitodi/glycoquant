# 라벨 N-글리칸 자동 분석기

Thermo `.raw` 또는 `.mzML`에서 **환원말단 표지(ProA 등) N-글리칸을 자동으로
동정·정량**해 엑셀로 내보내는 CLI 도구. **Xcalibur(윈도우 전용) 없이 macOS/Windows에서**
동작한다.

Xcalibur에서 스캔 값을 손으로 옮기고, 이론질량을 계산기 시트에 넣고, adduct 강도를
하나씩 찾아 합산하고, **반복 실험을 평균 ± SD로 취합**하던 수작업을 자동화한다.
파일 하나는 물론 **폴더째(반복 여러 개)** 넣으면 개별 결과와 취합 결과를 한 번에 낸다.

자세한 설치·실행·옵션·YAML 규칙 설명은 [USAGE.md](USAGE.md)를 참고한다.

> 방법론 근거: *Moon et al., "Structural, quantitative, and functional characterization
> of N-glycans in porcine pancreatin extract, lipase, and α-amylase",
> Enzyme and Microbial Technology 196 (2026) 110830* (LC-ESI-HCD-MS/MS, ProA 표지).

기본 설정은 **ProA N-글리칸**이지만 **`configs/*.yaml` 한 파일만 바꾸면**
다른 라벨(2AB, 2AA, free), 단당, adduct, 탐색범위로 재사용할 수 있다.

---

## 주요 기능

| 수작업 | 자동화 기능 | 확인 내용 |
|---|---|---|
| Xcalibur 스캔값 전사 | `Screening` 시트 자동 생성 | 204 oxonium 12/12 정확 일치, 전체 MS2 덤프 지원 |
| 이론 monoisotopic 질량 계산 | 설정 기반 질량 엔진 | 계산기 44개 [M+H]+ < 0.001 ppm |
| adduct m/z 나열 | H/Na/K × 1~3가 자동 생성 | 19종 일치 |
| adduct 강도 합산 | XIC 면적 자동 합산 | 논문식 면적 정량 |
| 상대 %, PPM 매칭, 유형 분류 | 자동 계산 | 유형분포 수작업과 일치 |
| 반복 평균 ± SD 취합 | 폴더 배치 `aggregated.xlsx` | 결측 0 처리, 검출빈도 표시 |
| 진단 규칙 기반 스크리닝 | `--targets` YAML 규칙 | `Screening` + `구조찾기` 시트 |

---

## 파이프라인

```text
.raw -> mzML -> MS1/MS2 파싱
  -> 설정의 단당·라벨로 조성 후보 생성
  -> H/Na/K adduct 이론 m/z 생성
  -> MS2 precursor ppm 매칭
  -> 진단 oxonium 이온 확인
  -> MS1 정밀질량 게이트
  -> XIC 면적 합산과 상대량(%)
  -> High-mannose / Hybrid / Complex 분류 + Oxford 명명
  -> 엑셀 리포트
```

---

## 실데이터 검증

**알츠하이머 두정엽 수작업 시트 대조**

- 공통 글리칸 40개, Pearson r = 0.94, 중앙오차 0.27 %p, Spearman 0.75
- 상위 글리칸 근접: FA3 21.6 %(수작업 19.8), M5 20.3 %(18.7), FM3 3.61 %(3.53)

**Porcine pancreatin 논문 PPE Table 1 대조**

- recall 24/32, FA2 17.4 %(논문 17.7), 유형분포 HM 30 / Hyb 3 / Cplx 68 %(논문 31/4/64)
- `--ms1-first` 사용 시 Neu5Gc 1.0 % -> 15.7 % 회수(논문 14.0)

**단위 테스트**

| 항목 | 기준 | 현재 |
|---|---|---|
| 질량엔진 (계산기 44개 [M+H]+) | < 0.001 ppm | 0.0007 ppm |
| 논문 진단이온 6종 | < 2 ppm | < 1 ppm |
| Oxford 명명 (논문 Table 1) | 일치 | 8/8 |
| 정량 순위 상관 (시트 조성, 현재 분석 경로) | Spearman >= 0.70 | +0.76 |

---

## 결과 엑셀

기본 분석은 다음 시트를 생성한다.

- `Glycans`: 동정·정량된 글리칸 목록(Oxford명, 조성, 유형, m/z, RT, Evidence, 상대%)
- `Summary`: 유형별(High-mannose/Hybrid/Complex) 및 시알/푸코 집계
- `Adducts`: 글리칸별 검출 adduct 상세
- `Screening`: MS2 스캔별 진단이온 관측 m/z, precursor, monoisotope, isolation, charge

`--targets` 진단 스크리닝은 별도 엑셀에 다음 시트를 생성한다.

- `Screening`: 채택 스캔, 코어이온 관측 m/z, ppm 오차, precursor, monoisotope, charge, Features
- `구조찾기`: monoisotope 기준 그룹핑, `SC no.` 역참조, adduct 수기 입력 열

---

## 구조

```text
glycan_analyze.py      # CLI 진입점
setup_thermo.py        # OS별 ThermoRawFileParser 설치
run.command / run.bat  # 비개발자용 더블클릭 래퍼
configs/               # 분석/진단 YAML 설정
glyco/                 # 라이브러리
  config.py            # YAML -> 설정 객체
  chem.py              # 질량/이온/진단이온
  compositions.py      # 조성 후보 생성
  identify.py          # precursor 매칭 + 진단 확인
  quantify.py          # XIC apex/area 정량
  aggregate.py         # 반복 취합
  targets.py           # --targets 진단 규칙 로딩/검증
  report.py            # 엑셀 출력
tests/                 # 회귀 테스트
```

---

## 알아둘 점

- 동정은 결정론적 m/z 매칭과 진단이온 확인을 기반으로 한다.
- 정량(%)은 자동 EIC 면적이므로 사람이 읽은 값과 절대값이 다를 수 있다. 순위와 패턴 중심으로 해석한다.
- 구조 분류와 Oxford 명명은 조성 기반 근사다. 이성질체와 안테나 위치는 구분하지 않는다.
- 절대정량(pmol)은 범위 밖이다. 표준품 검량선 기반 UPLC 형광 정량이 필요하다.

---

## 인용 / 라이선스

- 방법론 인용: Moon et al., *Enzyme and Microbial Technology* 196 (2026) 110830.
- 번들: [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) (자체 라이선스).
- 라이선스: [MIT](LICENSE)
