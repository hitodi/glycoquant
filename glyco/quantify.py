"""
정량 (XIC)
----------
타깃 m/z 집합에 대해 MS1 전체를 훑어 XIC(추출 이온 크로마토그램)를 만들고,
- apex : 최대 피크 높이
- area : apex RT 주변 구간의 면적(사다리꼴 적분)  <- 논문 §2.5 권장 방식
둘 다 제공한다. 면적적분은 (동정 후) 선택된 소수 타깃에만 적용한다.
"""

import numpy as np

# numpy 2.0 에서 trapz -> trapezoid 로 개명됨(양쪽 지원)
_trapz = getattr(np, "trapezoid", None) or np.trapz


def build_targets(candidates, adduct_ions, ion_mz_fn):
    """candidates × adducts -> (정렬된 mz ndarray, meta list)."""
    rows = []
    for ci, cand in enumerate(candidates):
        for adduct, z in adduct_ions:
            mz = ion_mz_fn(cand["neutral"], adduct, z)
            rows.append((mz, ci, adduct, z))
    rows.sort(key=lambda r: r[0])
    mz_sorted = np.array([r[0] for r in rows], dtype=np.float64)
    meta = [{"cand_idx": r[1], "adduct": r[2], "z": r[3], "mz": r[0]} for r in rows]
    return mz_sorted, meta


def extract_xic(msdata, mz_targets, ppm_tol=10.0, log=print):
    """
    선택된 타깃들의 XIC 강도행렬을 만든다.
    반환 dict:
      rt    : (n_scan,)          MS1 RT 축
      inten : (n_target, n_scan) 각 타깃의 스캔별 강도(매칭 없으면 0)
      mz    : (n_target, n_scan) 매칭된 실제 m/z (없으면 nan)
    """
    n_scan = len(msdata.ms1)
    nt = len(mz_targets)
    rt = np.array([s[0] for s in msdata.ms1], dtype=np.float64)
    inten = np.zeros((nt, n_scan), dtype=np.float64)
    mzmat = np.full((nt, n_scan), np.nan, dtype=np.float64)
    tol = mz_targets * ppm_tol * 1e-6

    for k, (_, mz, it) in enumerate(msdata.ms1):
        if mz.size == 0:
            continue
        idx = np.clip(np.searchsorted(mz, mz_targets), 1, mz.size - 1)
        left = idx - 1
        dl = np.abs(mz[left] - mz_targets)
        dr = np.abs(mz[idx] - mz_targets)
        use_left = dl <= dr
        nearest = np.where(use_left, left, idx)
        ndist = np.where(use_left, dl, dr)
        hit = ndist <= tol
        inten[hit, k] = it[nearest][hit]
        mzmat[hit, k] = mz[nearest][hit]
        if log and k % 500 == 0 and k:
            log(f"[정량] MS1 {k}/{n_scan} 스캔...")
    return {"rt": rt, "inten": inten, "mz": mzmat}


def apex(xic):
    """타깃별 apex 강도/ RT / 정밀 m/z."""
    inten, rt, mzmat = xic["inten"], xic["rt"], xic["mz"]
    ai = np.argmax(inten, axis=1)
    rows = np.arange(inten.shape[0])
    return inten[rows, ai], rt[ai], mzmat[rows, ai]


def area(xic, rt_window=0.5):
    """
    타깃별 면적: apex RT ± rt_window 구간을 사다리꼴 적분.
    apex 정밀 m/z 도 함께 반환(MS1 게이트용).
    """
    inten, rt, mzmat = xic["inten"], xic["rt"], xic["mz"]
    nt = inten.shape[0]
    ai = np.argmax(inten, axis=1)
    rows = np.arange(nt)
    apex_rt = rt[ai]
    apex_mz = mzmat[rows, ai]
    areas = np.zeros(nt, dtype=np.float64)
    for i in range(nt):
        if inten[i, ai[i]] <= 0:
            continue
        lo, hi = apex_rt[i] - rt_window, apex_rt[i] + rt_window
        sel = (rt >= lo) & (rt <= hi)
        if sel.sum() >= 2:
            areas[i] = _trapz(inten[i, sel], rt[sel])
        else:
            areas[i] = inten[i, ai[i]]   # 점 하나뿐이면 높이로 대체
    return areas, apex_rt, apex_mz


# --- 하위호환: 단일 함수형 apex 추출 ---
def xic_apex(msdata, mz_targets, ppm_tol=10.0, log=print):
    xic = extract_xic(msdata, mz_targets, ppm_tol=ppm_tol, log=log)
    return apex(xic)
