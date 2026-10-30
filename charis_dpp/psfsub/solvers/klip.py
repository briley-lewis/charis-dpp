"""KLIP (Soummer et al. 2012) with the Pueyo (2016) forward model.

Fixes vs. ``charis_adiklip.pro``: outlier pixels are NaN in the residual
instead of being left unsubtracted (#145); eigenvectors and eigenvalues are
ordered consistently for ``npca = 1`` (#146); degenerate eigenvalues are
skipped in the indirect term (#147); zone bookkeeping cannot desynchronise
(#144).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from charis_dpp.psfsub.common import outlier_mask

__all__ = ["klip_basis", "klip_zone"]


def klip_basis(refs: np.ndarray, npca: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """KL modes of the reference set ``refs`` ``(npix, nref)``.

    Returns ``(Z, lam, V)``: ``Z`` ``(npca, npix)`` orthonormal modes, the
    eigenvalues ``lam`` (descending) and eigenvectors ``V`` ``(nref, npca)``
    of the ``nref×nref`` covariance.
    """
    nref = refs.shape[1]
    cov = refs.T @ refs
    npca = max(1, min(npca, nref))
    lam, V = eigh(cov, subset_by_index=[nref - npca, nref - 1])
    lam, V = lam[::-1], V[:, ::-1]
    lam = np.where(lam > 0, lam, np.nan)
    Z = (V.T @ refs.T) / np.sqrt(lam)[:, None]
    return Z, lam, V


def klip_zone(
    sub: np.ndarray,
    *,
    ref_for: callable,
    npca: int,
    meansub: bool,
    zero: bool,
    sub_model: np.ndarray | None = None,
    min_refs: int = 2,
    min_pix: int = 5,
) -> np.ndarray:
    """KLIP-subtract every frame of one zone (``sub`` is ``(npix, nf)``; opt region == sub region).

    With ``sub_model`` the forward-modelled residual of the model (Pueyo 2016,
    over-subtraction + direct and indirect self-subtraction terms) is returned.
    """
    npix, nfiles = sub.shape
    fwd = sub_model is not None
    data = sub.copy()
    model = sub_model.copy() if fwd else None
    if meansub:
        data = data - np.nanmean(data, axis=0, keepdims=True)
        if fwd:
            model = model - np.nanmean(model, axis=0, keepdims=True)
    good = outlier_mask(data)
    out = np.full(sub.shape, np.nan)
    if good.sum() < min_pix:
        return out
    D = data[good]
    M = model[good] if fwd else None
    for nf in range(nfiles):
        indim = np.asarray(ref_for(nf), dtype=int)
        indim = indim[indim != nf]
        if indim.size < min_refs or np.isfinite(D[:, nf]).sum() < min_pix:
            continue
        R = D[:, indim]
        target = np.nan_to_num(D[:, nf])
        try:
            Z, lam, V = klip_basis(R, npca)
        except np.linalg.LinAlgError:
            continue
        ok = np.isfinite(lam)
        Z, lam, V = Z[ok], lam[ok], V[:, ok]
        if Z.shape[0] == 0:
            continue
        if not fwd:
            coeff = Z @ target
            ref = coeff @ Z
            diff = np.full(npix, np.nan)
            diff[good] = D[:, nf] - ref
        else:
            S = np.nan_to_num(M[:, nf])
            Rm = M[:, indim]
            c1 = indim.size
            cov = R.T @ R
            lam_full, V_full = eigh(cov)
            lam_full, V_full = lam_full[::-1], V_full[:, ::-1]
            pc_full = V_full.T @ R.T  # (c1, npix), unnormalised
            CdS = Rm.T @ R + R.T @ Rm  # (c1, c1)
            dZ = (V.T @ Rm.T) / np.sqrt(lam)[:, None]  # direct term
            for k in range(Z.shape[0]):
                vk = V[:, k]
                first = (-0.5 / lam[k]) * (vk @ CdS @ vk) * Z[k]
                second = np.zeros(Z.shape[1])
                for j in range(c1):
                    if (
                        j == k
                        or not (lam_full[j] > 0)
                        or abs(lam[k] - lam_full[j]) < 1e-12 * lam[k]
                    ):
                        continue
                    con = np.sqrt(lam_full[j] / lam[k]) / (lam[k] - lam_full[j])
                    second += (
                        con
                        * (V_full[:, j] @ CdS @ V_full[:, k])
                        * pc_full[j]
                        / np.sqrt(lam_full[j])
                    )
                dZ[k] = dZ[k] + first + second
            over = (Z @ S) @ Z
            self1 = (dZ @ target) @ Z
            self2 = (Z @ target) @ dZ
            diff = np.full(npix, np.nan)
            diff[good] = S - (over + self1 + self2)
        if zero:
            diff = diff - np.nanmedian(diff)
        out[:, nf] = diff
    return out
