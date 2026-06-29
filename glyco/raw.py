"""
.raw -> mzML 변환
-----------------
Thermo .raw 파일을 ThermoRawFileParser(번들)로 표준 mzML로 변환한다.
OS를 감지해 알맞은 실행 방식을 고른다.

  - Windows : vendor/thermo-win/ThermoRawFileParser.exe  (네이티브 실행)
  - macOS/Linux : vendor/thermo-<plat>/ThermoRawFileParser.dll  (dotnet 으로 실행)

번들이 없으면 setup_thermo.py 안내 메시지를 띄운다.
"""

import os
import platform
import shutil
import subprocess
import sys

VENDOR = os.path.join(os.path.dirname(__file__), "vendor")


def _plat_dir():
    sysname = platform.system()
    if sysname == "Windows":
        return "thermo-win", "exe"
    if sysname == "Darwin":
        arch = "arm64" if platform.machine() in ("arm64", "aarch64") else "x64"
        return f"thermo-osx-{arch}", "dll"
    return "thermo-linux", "dll"


def _find_parser():
    sub, kind = _plat_dir()
    base = os.path.join(VENDOR, sub)
    if kind == "exe":
        exe = os.path.join(base, "ThermoRawFileParser.exe")
        if os.path.isfile(exe):
            return ("native", exe)
    else:
        dll = os.path.join(base, "ThermoRawFileParser.dll")
        if os.path.isfile(dll):
            return ("dotnet", dll)
    return None


def _have_dotnet():
    return shutil.which("dotnet") is not None


def convert(raw_path: str, out_dir: str, indexed: bool = True, log=print) -> str:
    """
    raw_path -> mzML. 변환된 mzML 경로를 반환.
    """
    raw_path = os.path.abspath(raw_path)
    if not os.path.isfile(raw_path):
        raise FileNotFoundError(f"raw 파일을 찾을 수 없습니다: {raw_path}")
    os.makedirs(out_dir, exist_ok=True)

    found = _find_parser()
    if not found:
        sub, _ = _plat_dir()
        raise RuntimeError(
            f"ThermoRawFileParser 번들이 없습니다 (vendor/{sub}).\n"
            f"  먼저 'python setup_thermo.py' 를 실행해 변환기를 내려받으세요."
        )
    mode, target = found

    base = os.path.splitext(os.path.basename(raw_path))[0]
    out_mzml = os.path.join(out_dir, base + ".mzML")
    fmt = "2" if indexed else "1"  # 2=indexed mzML, 1=mzML

    if mode == "native":
        cmd = [target]
    else:
        if not _have_dotnet():
            raise RuntimeError(
                "이 OS에서는 .NET 런타임(dotnet)이 필요합니다.\n"
                "  https://dotnet.microsoft.com/download 에서 설치하세요."
            )
        cmd = ["dotnet", target]
    cmd += ["-i", raw_path, "-b", out_mzml, "-f", fmt, "-l", "3"]

    log(f"[변환] {os.path.basename(raw_path)} -> mzML ...")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.isfile(out_mzml):
        raise RuntimeError(
            "mzML 변환 실패:\n" + (proc.stderr or proc.stdout or "(no output)")
        )
    log(f"[변환] 완료 -> {out_mzml}")
    return out_mzml


def to_mzml(path: str, out_dir: str, log=print) -> str:
    """입력이 .raw 면 변환하고, 이미 .mzML 이면 그대로 사용."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mzml":
        return os.path.abspath(path)
    if ext == ".raw":
        return convert(path, out_dir, log=log)
    raise ValueError(f"지원하지 않는 입력 형식입니다: {ext} (.raw 또는 .mzML)")
