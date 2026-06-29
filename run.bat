@echo off
REM Windows 더블클릭 실행용. .raw 파일을 이 .bat 위로 끌어다 놓거나, 더블클릭 후 경로 입력.
cd /d "%~dp0"

where py >nul 2>nul && (set PY=py) || (set PY=python)

%PY% -c "import pyteomics, psims, numpy, openpyxl, lxml" 2>nul
if errorlevel 1 (
  echo [설치] 필요한 패키지를 설치합니다...
  %PY% -m pip install -r requirements.txt
)

REM 변환기(win 빌드)가 없으면 자동 설치
if not exist "glyco\vendor\thermo-win\ThermoRawFileParser.exe" (
  echo [설치] Thermo 변환기를 내려받습니다...
  %PY% setup_thermo.py
)

set RAW=%~1
if "%RAW%"=="" set /p RAW=".raw 파일 경로를 입력(또는 드래그)하세요: "

%PY% glycan_analyze.py %RAW%
echo.
pause
