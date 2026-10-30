"""LOCI least-squares reference construction (Lafrenière et al. 2007) with
SVD truncation and correlation-based frame selection (A-LOCI, Currie et al.
2012), and the perturbed-coefficient forward model (Currie et al. 2018).

Fixes vs. ``charis_adialoci.pro``: the outlier mask is applied on every path
(#132); the forward-model perturbation uses the same pixel set as the
covariance (#130/#131); ``nref = 0`` is handled (#133).
"""

from __future__ import annotations

import numpy as np

from charis_dpp.psfsub.common import outlier_mask
from charis_dpp.psfsub.frame_selection import most_correlated

__all__ = ["loci_solve", "loci_zone"]


def loci_solve(cov: np.ndarray, indim: np.ndarray, nf: int, svd_eps: float | None) -> np.ndarray:
    """Coefficients ``c`` minimising ``|target - Σ c_k ref_k|²`` from the covariance matrix."""
    A = cov[np.ix_(indim, indim)]
    b = cov[indim, nf]
    if svd_eps is not None and A.shape[0] > 1:
        return np.linalg.pinv(A, rcond=svd_eps) @ b
    return np.linalg.solve(A, b)


def loci_zone(
    opt: np.ndarray,
    sub: np.ndarray,
    *,
    ref_for: callable,
    nref: int,
    svd_eps: float | None,
    zero: bool,
    opt_model: np.ndarray | None = None,
    sub_model: np.ndarray | None = None,
    min_refs: int = 1,
    min_pix: int = 2,
) -> np.ndarray:
    """Subtract a LOCI reference from every frame in one zone.

    ``opt``/``sub`` are ``(npix, nf)`` matrices of the optimisation and
    subtraction regions; ``ref_for(nf)`` returns the candidate reference
    frame indices for frame ``nf``.  With ``opt_model``/``sub_model`` the
    forward-modelled residual of the model is returned instead.  Returns
    ``(npix_sub, nf)`` (NaN where no reference could be built).
    """
    npix, nfiles = sub.shape
    fwd = opt_model is not None
    good_opt = outlier_mask(opt)
    if good_opt.sum() < 5:
        return np.full(sub.shape, np.nan)
    X = opt[good_opt]
    Xm = opt_model[good_opt] if fwd else None
    cov = X.T @ X  # (nf, nf) data covariance
    out = np.full(sub.shape, np.nan)
    for nf in range(nfiles):
        indim = np.asarray(ref_for(nf), dtype=int)
        indim = indim[indim != nf]
        if indim.size < min_refs:
            continue
        indim = most_correlated(cov, nf, indim, nref)
        good_sub = np.isfinite(sub[:, nf])
        if good_sub.sum() < min_pix:
            continue
        try:
            c = loci_solve(cov, indim, nf, svd_eps)
        except np.linalg.LinAlgError:
            continue
        if not fwd:
            ref = np.nansum(sub[:, indim] * c[None, :], axis=1)
            diff = sub[:, nf] - ref
        else:
            # perturbation of the coefficients by the model (first-order term)
            refpert = Xm[:, indim] @ c
            pert2 = Xm[:, nf] - refpert
            b = X[:, indim].T @ pert2
            try:
                dc = loci_solve_matrix(cov, indim, b, svd_eps)
            except np.linalg.LinAlgError:
                continue
            dc = np.where(np.isfinite(dc), dc, 0.0)
            ref = (sub_model[:, indim] @ c) + (sub[:, indim] @ dc) + (sub_model[:, indim] @ dc)
            diff = sub_model[:, nf] - ref
        diff = np.where(good_sub, diff, np.nan)
        if zero:
            diff = diff - np.nanmedian(diff)
        out[:, nf] = diff
    return out


def loci_solve_matrix(
    cov: np.ndarray, indim: np.ndarray, b: np.ndarray, svd_eps: float | None
) -> np.ndarray:
    A = cov[np.ix_(indim, indim)]
    if svd_eps is not None and A.shape[0] > 1:
        return np.linalg.pinv(A, rcond=svd_eps) @ b
    return np.linalg.solve(A, b)
