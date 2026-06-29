"""
정량 (XIC 기반)
---------------
타깃 m/z 목록에 대해 MS1 전체를 훑어 각 타깃의 XIC 피크 높이(apex)와
그 RT를 구한다. (정량 방식은 기존 시트와 동일: adduct별 강도를 합산)

성능: MS1 스펙트럼마다 벡터화된 searchsorted로 모든 피크를 한 번에 매칭하므로
      수만 스캔도 수 초~수십 초에 처리된다.
"""

import numpy as np


def build_targets(candidates, adduct_ions, ion_mz_fn):
    """
    candidates: compositions.generate() 결과 리스트
    adduct_ions: [(adduct_tuple, z), ...]
    ion_mz_fn: (neutral, adduct, z) -> m/z
    반환: (mz_sorted ndarray, meta list) — meta[i] = dict(cand_idx, adduct, z, mz)
    """
    rows = []
    for ci, cand in enumerate(candidates):
        for adduct, z in adduct_ions:
            mz = ion_mz_fn(cand["neutral"], adduct, z)
            rows.append((mz, ci, adduct, z))
    rows.sort(key=lambda r: r[0])
    mz_sorted = np.array([r[0] for r in rows], dtype=np.float64)
    meta = [{"cand_idx": r[1], "adduct": r[2], "z": r[3], "mz": r[0]} for r in rows]
    return mz_sorted, meta


def xic_apex(msdata, mz_targets, ppm_tol=10.0, log=print):
    """
    각 타깃 m/z의 XIC 최대강도(apex height)와 그 RT, 그리고 apex 시점에
    실제로 매칭된 MS1 피크의 m/z(정밀질량)를 구한다.
    반환: (apex_intensity, apex_rt, apex_mz) — 모두 mz_targets 와 같은 길이
    """
    nt = len(mz_targets)
    best_int = np.zeros(nt, dtype=np.float64)
    best_rt = np.full(nt, np.nan, dtype=np.float64)
    best_mz = np.full(nt, np.nan, dtype=np.float64)
    tol = mz_targets * ppm_tol * 1e-6   # 타깃별 절대 허용오차

    for k, (rt, mz, inten) in enumerate(msdata.ms1):
        if mz.size == 0:
            continue
        # 각 타깃에 대해 스펙트럼에서 가장 가까운 피크를 찾는다
        idx = np.searchsorted(mz, mz_targets)
        idx = np.clip(idx, 1, mz.size - 1)
        left = idx - 1
        # 왼/오 후보 중 더 가까운 쪽 선택
        dl = np.abs(mz[left] - mz_targets)
        dr = np.abs(mz[idx] - mz_targets)
        use_left = dl <= dr
        nearest = np.where(use_left, left, idx)
        ndist = np.where(use_left, dl, dr)
        peak_int = inten[nearest]
        hit = (ndist <= tol)
        # apex 갱신
        upd = hit & (peak_int > best_int)
        best_int[upd] = peak_int[upd]
        best_rt[upd] = rt
        best_mz[upd] = mz[nearest][upd]
        if log and k % 500 == 0 and k:
            log(f"[정량] MS1 {k}/{len(msdata.ms1)} 스캔...")
    return best_int, best_rt, best_mz
