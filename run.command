#!/bin/bash
# macOS 더블클릭 실행용. .raw 파일을 이 아이콘 위로 끌어다 놓거나, 더블클릭 후 경로 입력.
cd "$(dirname "$0")"

PY=python3
command -v "$PY" >/dev/null 2>&1 || PY=python

# 최초 1회 의존성 설치 확인
"$PY" -c "import pyteomics, psims, numpy, openpyxl, lxml" 2>/dev/null || {
  echo "[설치] 필요한 패키지를 설치합니다..."
  "$PY" -m pip install -r requirements.txt
}

RAW="$1"
if [ -z "$RAW" ]; then
  read -r -p ".raw 파일 경로를 입력(또는 드래그)하세요: " RAW
  RAW="${RAW%\"}"; RAW="${RAW#\"}"   # 따옴표 제거
fi

"$PY" glycan_analyze.py "$RAW"
echo
read -r -p "끝났습니다. Enter 를 누르면 닫힙니다."
