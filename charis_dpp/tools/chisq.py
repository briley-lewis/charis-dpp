"""χ² of a model spectrum against data with an analytic scale factor
(``subroutines/chisq_calcs/*.pro``).

The IDL routines brute-forced the scale factor α over 10⁴–10⁵ steps; here
the weighted least-squares α is solved analytically (or by a 1-D bounded
minimisation when model errors scale with α).  With a covariance matrix,
excluded channels are dropped from the matrix (the IDL version zeroed the
residual but kept the full matrix, #295).  ``ndof`` counts the photometric
points and the fitted parameter (#296).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

__all__ = ["ChiSq", "chisq", "chisq_covar", "binned_chisq"]


@dataclass
class ChiSq:
    chisq: float
    ndof: int
    alpha: float

    @property
    def reduced(self) -> float:
        return self.chisq / max(self.ndof, 1)


def _minimise_log(f, a0: float, bounds) -> tuple[float, float]:
    """1-D bounded minimisation of ``f(alpha)`` in log(alpha): the scale factor can be ~1e-20
    (model surface flux vs. observed mJy), far below the absolute tolerance of a linear search."""
    a0 = abs(a0) if a0 != 0 and np.isfinite(a0) else 1.0
    if bounds is None:
        lo, hi = np.log10(a0) - 4.0, np.log10(a0) + 4.0
    else:
        lo, hi = np.log10(min(bounds)), np.log10(max(bounds))
    r = minimize_scalar(
        lambda la: f(10.0**la), bounds=(lo, hi), method="bounded", options={"xatol": 1e-7}
    )
    return float(10.0**r.x), float(r.fun)


def _prep(obs, eobs, model, exclude):
    obs = np.asarray(obs, dtype=float)
    eobs = np.asarray(eobs, dtype=float)
    model = np.asarray(model, dtype=float)
    good = np.isfinite(obs) & np.isfinite(eobs) & np.isfinite(model) & (eobs > 0)
    if exclude is not None:
        good[np.asarray(exclude, dtype=int)] = False
    return obs, eobs, model, good


def chisq(
    obs,
    eobs,
    model,
    *,
    err_model=None,
    exclude=None,
    phot=None,
    ephot=None,
    modelphot=None,
    alpha_bounds=None,
) -> ChiSq:
    """Diagonal χ² with the best scale factor α applied to the model."""
    obs, eobs, model, good = _prep(obs, eobs, model, exclude)
    o, e, m = obs[good], eobs[good], model[good]
    pm = np.asarray(phot, float) if phot is not None else np.zeros(0)
    pe = np.asarray(ephot, float) if phot is not None else np.zeros(0)
    pmod = np.asarray(modelphot, float) if phot is not None else np.zeros(0)
    if err_model is None:
        w = 1 / np.concatenate([e, pe]) ** 2
        y = np.concatenate([o, pm])
        x = np.concatenate([m, pmod])
        alpha = float(np.sum(w * x * y) / np.sum(w * x * x))
        chi = float(np.sum(w * (y - alpha * x) ** 2))
    else:
        em = np.asarray(err_model, float)[good]

        def f(a):
            v = e**2 + a**2 * em**2
            c = np.sum((o - a * m) ** 2 / v)
            if pm.size:
                c += np.sum(((pm - a * pmod) / pe) ** 2)
            return c

        a0 = float(np.sum(o * m / e**2) / np.sum(m * m / e**2))  # analytic start
        alpha, chi = _minimise_log(f, a0, alpha_bounds)
    return ChiSq(chi, int(o.size + pm.size - 1), alpha)


def chisq_covar(
    obs,
    eobs,
    model,
    speccovar,
    *,
    err_model=None,
    exclude=None,
    phot=None,
    ephot=None,
    modelphot=None,
) -> ChiSq:
    """χ² with a full spectral covariance (correlation ψ scaled by the errors)."""
    obs, eobs, model, good = _prep(obs, eobs, model, exclude)
    idx = np.flatnonzero(good)
    psi = np.asarray(speccovar, float)[np.ix_(idx, idx)]
    o, e, m = obs[idx], eobs[idx], model[idx]
    em = np.asarray(err_model, float)[idx] if err_model is not None else None
    pm = np.asarray(phot, float) if phot is not None else np.zeros(0)
    pe = np.asarray(ephot, float) if phot is not None else np.zeros(0)
    pmod = np.asarray(modelphot, float) if phot is not None else np.zeros(0)

    def f(a):
        cov = psi * np.outer(e, e)
        if em is not None:
            cov = cov + np.diag(a**2 * em**2)
        d = o - a * m
        try:
            c = float(d @ np.linalg.solve(cov, d))
        except np.linalg.LinAlgError:
            c = float(d @ np.linalg.pinv(cov) @ d)
        if pm.size:
            c += float(np.sum(((pm - a * pmod) / pe) ** 2))
        return c

    a0 = float(np.sum(o * m / e**2) / np.sum(m * m / e**2))
    alpha, chi = _minimise_log(f, a0, None)
    return ChiSq(chi, int(o.size + pm.size - 1), alpha)


def binned_chisq(chisq_map: np.ndarray, fwhm: float) -> tuple[float, int, np.ndarray]:
    """Reduced χ² of a per-pixel χ² map binned to resolution elements (block means, not nearest-neighbour sampling)."""
    m = np.asarray(chisq_map, float)
    b = max(int(round(fwhm)), 1)
    ny, nx = (m.shape[0] // b) * b, (m.shape[1] // b) * b
    blocks = m[:ny, :nx].reshape(ny // b, b, nx // b, b)
    binned = np.nanmean(blocks, axis=(1, 3))
    good = np.isfinite(binned) & (binned != 0)
    n = int(good.sum())
    total = float(np.nansum(binned[good]))
    return total / max(n - 1, 1), n, binned
