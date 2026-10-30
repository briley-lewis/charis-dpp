"""Aperture photometry (astrolib ``aper`` with ``/EXACT``, ``charis_myaper``, ``sumaper_im``)."""

from __future__ import annotations

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAnnulus, CircularAperture
from scipy.signal import fftconvolve

__all__ = ["aperture_weights", "aper_idl", "aperture_sum", "sumaper_im"]


def aperture_weights(shape, x: float, y: float, r: float) -> np.ndarray:
    """Exact fractional pixel coverage of a circle of radius ``r`` at ``(x, y)`` (full-frame array)."""
    ap = CircularAperture((x, y), r=max(r, 1e-6))
    m = ap.to_mask(method="exact")
    return m.to_image(shape[-2:])


def aperture_sum(im: np.ndarray, x: float, y: float, r: float, *, nan: bool = True) -> float:
    """Sum of ``im`` inside the exact aperture; NaN pixels count as 0 when ``nan``."""
    w = aperture_weights(im.shape, x, y, r)
    data = np.nan_to_num(im) if nan else im
    return float(np.sum(w * data))


def aper_idl(
    im: np.ndarray,
    x: float,
    y: float,
    r: float,
    skyrad: tuple[float, float] | None = None,
    *,
    setskyval: float | None = 0.0,
    phpadu: float = 1.0,
    sky_method: str = "median",
) -> tuple[float, float, float, float]:
    """``aper, im, x, y, flux, eflux, sky, skyerr, r, skyrad, /exact, /nan, /flux [, setskyval=]``.

    Returns ``(flux, eflux, sky, skyerr)``.  ``setskyval=None`` estimates the
    sky in the annulus ``skyrad=(rin, rout)`` with ``sky_method`` ``"median"``
    (plain median) or ``"mmm"`` (sigma-clipped, close to astrolib ``mmm``).
    """
    im = np.asarray(im, dtype=float)
    w = aperture_weights(im.shape, x, y, r)
    good = np.isfinite(im)
    area = float(np.sum(w))
    if setskyval is not None:
        sky, skyerr, nsky = float(setskyval), 0.0, 1
    else:
        if skyrad is None:
            raise ValueError("aper_idl needs skyrad when setskyval is None")
        ann = CircularAnnulus((x, y), r_in=skyrad[0], r_out=skyrad[1])
        am = ann.to_mask(method="center").to_image(im.shape[-2:]) > 0
        vals = im[am & good]
        if vals.size < 2:
            sky, skyerr, nsky = 0.0, 0.0, 1
        elif sky_method == "mmm":
            _, sky, skyerr = sigma_clipped_stats(vals, sigma=3.0)
            nsky = vals.size
        else:
            sky, skyerr, nsky = float(np.median(vals)), float(np.std(vals)), vals.size
    flux = float(np.sum(w * np.where(good, im - sky, 0.0)))
    var = max(flux, 0.0) / phpadu + area * skyerr**2 + area**2 * skyerr**2 / max(nsky, 1)
    return flux, float(np.sqrt(var)), float(sky), float(skyerr)


def sumaper_im(
    im: np.ndarray, r: float, rmax: float | None = None, center=None, *, nan: bool = True
) -> np.ndarray:
    """Image whose value at each pixel is the sum of ``im`` in a radius-``r`` aperture centred there.

    Pixels farther than ``rmax`` from ``center`` are 0, as in ``sumaper_im.pro``.
    NaN input pixels contribute 0.
    """
    from charis_dpp.geometry.grids import radial_grid

    im = np.asarray(im, dtype=float)
    k = int(np.ceil(r)) + 1
    kern = aperture_weights((2 * k + 1, 2 * k + 1), k, k, r)
    data = np.nan_to_num(im) if nan else im
    out = fftconvolve(data, kern, mode="same")
    if rmax is not None:
        out[radial_grid(im.shape, center) > rmax] = 0.0
    return out
