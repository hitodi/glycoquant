"""
mzML 파싱
---------
MS1 스펙트럼(정량용 XIC 재료)과 MS2 precursor(구조 확인용)를 읽어들인다.
대용량(수만 스캔)이라 MS1은 numpy 배열로 메모리에 적재한다(보통 수백 MB 이내).
"""

import numpy as np
from pyteomics import mzml


def _float_or_none(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class MsData:
    """파싱된 MS1/MS2 데이터를 담는 컨테이너."""

    def __init__(self):
        self.ms1 = []        # [(rt, mz_array, intensity_array), ...]  RT 오름차순
        self.ms2 = []        # [(rt, precursor_mz, charge_or_None, scan), ...]
        self.ms2_peaks = {}  # scan -> (mz_array, intensity_array)  (keep_ms2_peaks=True 일 때)
        self.ms2_info = {}   # scan -> precursor 부가정보(selected/isolation/charge intensity 등)

    @property
    def rt_range(self):
        if not self.ms1:
            return (0.0, 0.0)
        return (self.ms1[0][0], self.ms1[-1][0])


def parse(mzml_path: str, keep_ms2_peaks: bool = False, keep_ms1: bool = True, log=print) -> MsData:
    """
    keep_ms2_peaks=True 면 MS2 단편 피크도 보관(진단 oxonium 확인용, 메모리 추가).
    keep_ms1=False 면 스크리닝 전용으로 MS1 배열 적재를 건너뛴다.
    """
    data = MsData()
    n = 0
    for spec in mzml.read(mzml_path):
        n += 1
        level = spec.get("ms level")
        scan_info = spec["scanList"]["scan"][0]
        rt = float(scan_info["scan start time"])  # 분 단위(보통)
        if level == 1 and keep_ms1:
            mz = np.asarray(spec["m/z array"], dtype=np.float64)
            inten = np.asarray(spec["intensity array"], dtype=np.float64)
            data.ms1.append((rt, mz, inten))
        elif level == 2:
            try:
                prec = spec["precursorList"]["precursor"][0]
                ion = prec["selectedIonList"]["selectedIon"][0]
                pmz = float(ion.get("selected ion m/z"))
                z = ion.get("charge state")
                z = int(z) if z is not None else None
            except (KeyError, IndexError, TypeError):
                continue
            scan = spec.get("id", "").split("scan=")[-1]
            data.ms2.append((rt, pmz, z, scan))
            isolation = prec.get("isolationWindow", {})
            data.ms2_info[scan] = {
                "selected_mz": pmz,
                # Thermo 가 계산한 monoisotopic m/z(스캔 userParam) — 스크리닝 시트의 monoisotope 열
                "monoisotopic_mz": _float_or_none(scan_info.get("[Thermo Trailer Extra]Monoisotopic M/Z:")),
                "isolation_target_mz": _float_or_none(isolation.get("isolation window target m/z")),
                "isolation_lower_offset": _float_or_none(isolation.get("isolation window lower offset")),
                "isolation_upper_offset": _float_or_none(isolation.get("isolation window upper offset")),
                "charge": z,
                "peak_intensity": _float_or_none(ion.get("peak intensity")),
            }
            if keep_ms2_peaks:
                data.ms2_peaks[scan] = (
                    np.asarray(spec["m/z array"], dtype=np.float64),
                    np.asarray(spec["intensity array"], dtype=np.float64),
                )
        if log and n % 2000 == 0:
            log(f"[파싱] {n} 스캔...")
    data.ms1.sort(key=lambda x: x[0])
    data.ms2.sort(key=lambda x: x[0])
    log(f"[파싱] 완료 — MS1 {len(data.ms1)}개, MS2 {len(data.ms2)}개"
        + (f", MS2피크 {len(data.ms2_peaks)}개 보관" if keep_ms2_peaks else ""))
    return data
