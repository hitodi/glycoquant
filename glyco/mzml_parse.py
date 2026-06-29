"""
mzML 파싱
---------
MS1 스펙트럼(정량용 XIC 재료)과 MS2 precursor(구조 확인용)를 읽어들인다.
대용량(수만 스캔)이라 MS1은 numpy 배열로 메모리에 적재한다(보통 수백 MB 이내).
"""

import numpy as np
from pyteomics import mzml


class MsData:
    """파싱된 MS1/MS2 데이터를 담는 컨테이너."""

    def __init__(self):
        self.ms1 = []        # [(rt, mz_array, intensity_array), ...]  RT 오름차순
        self.ms2 = []        # [(rt, precursor_mz, charge_or_None, scan), ...]
        self.ms2_peaks = {}  # scan -> (mz_array, intensity_array)  (keep_ms2_peaks=True 일 때)

    @property
    def rt_range(self):
        if not self.ms1:
            return (0.0, 0.0)
        return (self.ms1[0][0], self.ms1[-1][0])


def parse(mzml_path: str, keep_ms2_peaks: bool = False, log=print) -> MsData:
    """
    keep_ms2_peaks=True 면 MS2 단편 피크도 보관(진단 oxonium 확인용, 메모리 추가).
    """
    data = MsData()
    n = 0
    for spec in mzml.read(mzml_path):
        n += 1
        level = spec.get("ms level")
        scan_info = spec["scanList"]["scan"][0]
        rt = float(scan_info["scan start time"])  # 분 단위(보통)
        if level == 1:
            mz = np.asarray(spec["m/z array"], dtype=np.float64)
            inten = np.asarray(spec["intensity array"], dtype=np.float64)
            data.ms1.append((rt, mz, inten))
        elif level == 2:
            try:
                ion = spec["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
                pmz = float(ion.get("selected ion m/z"))
                z = ion.get("charge state")
                z = int(z) if z is not None else None
            except (KeyError, IndexError, TypeError):
                continue
            scan = spec.get("id", "").split("scan=")[-1]
            data.ms2.append((rt, pmz, z, scan))
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
