"""Signal-to-noise maps (``charis_snratio_sub.pro``).

The image is convolved with a circular aperture of radius ``fwhm/2``; the
noise at each radius is the robust (MAD/0.6745) azimuthal scatter of that
aperture-summed image, optionally corrected for small-sample statistics
(Mawet et al. 2014, Student-t), and optionally with the azimuthal median
removed first (``zero``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

from charis_dpp.geometry.apertures import aper_idl, sumaper_im
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.centroid import gcntrd
from charis_dpp.image.filters import idl_filter_image, radial_profile_2d

__all__ = ["SNRResult", "snr_map", "small_sample_correction"]


@dataclass
class SNRResult:
    snrmap: np.ndarray
    noisemap: np.ndarray
    convimage: np.ndarray
    snr: float = np.nan  # at the requested position
    position: tuple[float, float] | None = None
    flux: float = np.nan


def small_sample_correction(
    shape, fwhm: float, center=None, *, with_source: bool = False, extra_masked=None
) -> np.ndarray:
    """Per-pixel factor multiplying the noise for the finite number of resolution elements at each radius."""
    r = radial_grid(shape, center)
    radpoint = np.arange(1, max(shape[-2:]) // 2 + 1, dtype=float)
    nres = 2 * np.pi * radpoint / fwhm
    lost = 2.0 if with_source else 1.0
    if extra_masked is not None:
        nres = nres - np.interp(radpoint, np.arange(1, len(extra_masked) + 1), extra_masked)
    dof = nres - lost
    corr = np.full(radpoint.size, 999.0)
    ok = dof > 0
    # t_cvf(2.867e-7, dof): the 5-sigma Gaussian tail probability as a t-cutoff
    corr[ok] = stats.t.isf(2.867e-7, dof[ok]) / 5.0 * np.sqrt(1 + 1 / dof[ok])
    return np.interp(r, radpoint, corr, left=999.0, right=corr[-1])


def snr_map(
    image: np.ndarray,
    fwhm: float,
    rmax: float,
    *,
    coord: tuple[float, float] | None = None,
    center: tuple[float, float] | None = None,
    exc: float = 1.5,
    highpass: bool = False,
    highpass_factor: float = 1.0,
    fixpos: bool = False,
    zero: bool = False,
    finite: bool = False,
    roi_exclude: np.ndarray | None = None,
) -> SNRResult:
    """Port of ``charis_snratio_sub``; see module docstring.

    ``coord`` gives a source position: the noise is then measured with the
    source masked (radius ``exc·fwhm/2``) and ``result.snr`` is the aperture
    flux at the (re-centroided, unless ``fixpos``) position divided by the
    local noise.  ``roi_exclude`` masks pixels (e.g. a disk) from the noise.
    """
    im = np.asarray(image, dtype=float)
    if highpass:
        im = im - idl_filter_image(im, median=int(round(5 * fwhm * highpass_factor)))
    xc, yc = image_center(im.shape) if center is None else center
    immask = np.ones_like(im)
    if roi_exclude is not None:
        immask[roi_exclude] = np.nan
    corr = (
        small_sample_correction(im.shape, fwhm, (xc, yc), with_source=coord is not None)
        if finite
        else 1.0
    )

    im1 = sumaper_im(im, fwhm / 2.0, rmax, (xc, yc))

    def noise_of(a: np.ndarray) -> np.ndarray:
        a = a.copy()
        if zero:
            a = a - radial_profile_2d(a * immask, 1.0, 0.0, rmax, (xc, yc))
        n = radial_profile_2d(np.abs(a * immask) / 0.6745, 1.0, 0.0, rmax, (xc, yc))
        n = n * corr
        return np.where(n > 0, n, np.nan)  # no azimuthal scatter → no SNR, not an infinite one

    noise = noise_of(im1)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = im1 / noise
    res = SNRResult(snrmap=snr, noisemap=noise, convimage=im1)
    if coord is None:
        return res
    xg, yg = coord
    if fixpos:
        xf, yf = float(xg), float(yg)
    else:
        xf, yf = gcntrd(im, xg, yg, fwhm)
        if not np.isfinite(xf):
            xf, yf = float(xg), float(yg)
    im2 = im.copy()
    im2[radial_grid(im.shape, (xf, yf)) <= exc * fwhm / 2.0] = np.nan
    im2 = sumaper_im(im2, fwhm / 2.0, rmax, (xc, yc))
    noise2 = noise_of(im2)
    with np.errstate(divide="ignore", invalid="ignore"):
        res.snrmap = im1 / noise2
    res.noisemap = noise2
    flux, *_ = aper_idl(im, xf, yf, fwhm * 0.5, setskyval=0.0)
    iy, ix = int(round(yf)), int(round(xf))
    if 0 <= iy < im.shape[0] and 0 <= ix < im.shape[1] and np.isfinite(flux) and noise2[iy, ix] > 0:
        res.snr = flux / noise2[iy, ix]
    res.flux = flux
    res.position = (xf, yf)
    return res
