"""
진단(oxonium) 이온 확인
-----------------------
MS2 단편 스펙트럼에 글리칸 공통 진단이온(예: 204.0867 HexNAc, 441.2708 ProA-HexNAc)이
있는지 확인한다. 화학적 노이즈 매칭을 걸러 '글리칸 vs 비글리칸'을 가린다.
(주의: oxonium 은 거의 모든 N-글리칸에 공통이라 '조성 vs 조성' 구분은 못 함.)
"""

import numpy as np


def present_ions(peaks_mz, diag_table, ppm_tol=20.0, min_rel_intensity=0.0, intensities=None):
    """
    peaks_mz : MS2 단편 m/z 배열(정렬 가정 X)
    diag_table : [(name, mz), ...]
    반환: 존재하는 진단이온 name 집합
    """
    if peaks_mz is None or len(peaks_mz) == 0:
        return set()
    mz = np.asarray(peaks_mz)
    order = np.argsort(mz)
    mz = mz[order]
    if intensities is not None and min_rel_intensity > 0:
        it = np.asarray(intensities)[order]
        thr = it.max() * min_rel_intensity
    else:
        it = None
        thr = 0
    found = set()
    for name, target in diag_table:
        idx = np.searchsorted(mz, target)
        for j in (idx, idx - 1):
            if 0 <= j < mz.size:
                if abs(mz[j] - target) / target * 1e6 <= ppm_tol:
                    if it is None or it[j] >= thr:
                        found.add(name)
                        break
    return found


def confirm_scan(msdata, scan, diag_table, required, ppm_tol=20.0):
    """해당 MS2 스캔이 required 진단이온을 모두 포함하면 True."""
    peaks = msdata.ms2_peaks.get(scan)
    if peaks is None:
        return None   # 피크 정보 없음(보관 안 함) -> 판단 불가
    found = present_ions(peaks[0], diag_table, ppm_tol=ppm_tol)
    return set(required).issubset(found)
