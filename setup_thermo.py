#!/usr/bin/env python3
"""
ThermoRawFileParser(변환기) 자동 설치
-------------------------------------
현재 OS에 맞는 ThermoRawFileParser 빌드를 내려받아 glyco/vendor/ 에 푼다.
macOS(arm64)용은 저장소에 이미 포함되어 있으므로 보통 윈도우 랩 PC에서만 필요.

사용:  python setup_thermo.py
"""

import io
import os
import platform
import sys
import urllib.request
import zipfile

RELEASE = "v.2.0.0-dev"
BASE = f"https://github.com/CompOmics/ThermoRawFileParser/releases/download/{RELEASE}"
VENDOR = os.path.join(os.path.dirname(__file__), "glyco", "vendor")

ASSETS = {
    "Windows": ("ThermoRawFileParser-{r}-win.zip", "thermo-win"),
    "Darwin-arm64": ("ThermoRawFileParser-{r}-osx-arm64.zip", "thermo-osx-arm64"),
    "Darwin-x86_64": ("ThermoRawFileParser-{r}-osx.zip", "thermo-osx-x64"),
    "Linux": ("ThermoRawFileParser-{r}-linux.zip", "thermo-linux"),
}


def key():
    s = platform.system()
    if s == "Darwin":
        return f"Darwin-{'arm64' if platform.machine() in ('arm64','aarch64') else 'x86_64'}"
    return s


def main():
    k = key()
    if k not in ASSETS:
        print(f"[오류] 지원하지 않는 OS: {k}")
        return 2
    asset, sub = ASSETS[k]
    dst = os.path.join(VENDOR, sub)
    if os.path.isdir(dst) and any(f.startswith("ThermoRawFileParser") for f in os.listdir(dst)):
        print(f"[설치] 이미 존재합니다: {dst}")
        return 0
    url = f"{BASE}/{asset.format(r=RELEASE)}"
    print(f"[다운로드] {url}")
    with urllib.request.urlopen(url) as resp:
        buf = resp.read()
    print(f"[압축해제] -> {dst}")
    os.makedirs(dst, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(buf)) as z:
        z.extractall(dst)
    # 네이티브 실행파일 권한
    for name in ("ThermoRawFileParser", "ThermoRawFileParser.exe"):
        p = os.path.join(dst, name)
        if os.path.isfile(p):
            os.chmod(p, 0o755)
    print("[완료] 변환기 설치 완료.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
