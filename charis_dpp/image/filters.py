"""Spatial filters and profiles.

* :func:`median_highpass` / :func:`idl_filter_image` — ``filter_image(median=N)``
  semantics (border restored to the input, as IDL does without ``/ALL_PIXELS``).
* :func:`radial_profile` / :func:`radial_profile_2d` — ``profrad_tc.pro``
  (the swapped ``interpol`` arguments at line 88 are fixed).
* :func:`azimuthal_profile` — ``profang_tc.pro`` honouring ``r1``/``r2``.
* :func:`fixpix` — ``fixpix_rs.pro`` with the sigma-buffer bug fixed.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.interpolate import CubicSpline

from charis_dpp.geometry.grids import image_center, radial_grid

__all__ = [
    "idl_filter_image",
    "median_highpass",
    "idl_smooth",
    "radial_profile",
    "radial_profile_2d",
    "azimuthal_profile",
    "fixpix",
]


def idl_filter_image(
    im: np.ndarray, median: int | None = None, smooth: int | None = None, all_pixels: bool = False
) -> np.ndarray:
    """Astrolib ``filter_image``: median (then boxcar) filter; NaNs ignored in the median.

    Without ``all_pixels`` the outer ``N//2`` pixels keep their input values,
    reproducing the default IDL behaviour the pipeline relied on.
    """
    im = np.asarray(im, dtype=float)
    out = im.copy()
    if median and median > 1:
        n = int(median)
        if np.isnan(im).any():
            out = _nan_median_filter(im, n)
        else:
            out = ndimage.median_filter(im, size=n, mode="nearest")
        if not all_pixels:
            b = n // 2
            out = _restore_border(out, im, b)
    if smooth and smooth > 1:
        n = int(smooth)
        sm = ndimage.uniform_filter(np.nan_to_num(out), size=n, mode="nearest")
        if not all_pixels:
            sm = _restore_border(sm, out, n // 2)
        out = sm
    return out


def _restore_border(filtered: np.ndarray, original: np.ndarray, b: int) -> np.ndarray:
    if b <= 0:
        return filtered
    out = filtered.copy()
    out[:b, :] = original[:b, :]
    out[-b:, :] = original[-b:, :]
    out[:, :b] = original[:, :b]
    out[:, -b:] = original[:, -b:]
    return out


def _nan_median_filter(im: np.ndarray, n: int) -> np.ndarray:
    """Median filter ignoring NaN (generic_filter is slow but exact)."""
    return ndimage.generic_filter(im, np.nanmedian, size=n, mode="nearest")


def median_highpass(im: np.ndarray, box: int) -> np.ndarray:
    """``im - filter_image(im, median=box)`` — the pipeline's spatial filter (medbox)."""
    return im - idl_filter_image(im, median=box)


def idl_smooth(im: np.ndarray, width: int) -> np.ndarray:
    """IDL ``smooth(im, width)`` (boxcar, edges untouched)."""
    w = int(width)
    if w <= 1:
        return np.asarray(im, dtype=float).copy()
    sm = ndimage.uniform_filter(np.asarray(im, dtype=float), size=w, mode="nearest")
    return _restore_border(sm, np.asarray(im, dtype=float), w // 2)


def radial_profile(
    im: np.ndarray,
    res: float = 1.0,
    rmin: float = 0.0,
    rmax: float | None = None,
    center: tuple[float, float] | None = None,
    *,
    statistic: str = "median",
) -> tuple[np.ndarray, np.ndarray]:
    """Azimuthal median (or max) in radial bins of width ``res`` (``profrad_tc`` ``p1d``).

    Returns ``(radii, profile)``; bins without finite pixels are NaN.
    """
    im = np.asarray(im, dtype=float)
    r = radial_grid(im.shape, center)
    if rmax is None:
        rmax = float(r.max())
    rr = np.round(r / res) * res if res > 0 else r
    sel = (r >= rmin) & (r <= rmax)
    radii = np.unique(rr[sel])
    bins = np.searchsorted(radii, rr[sel])
    vals = im[sel]
    prof = np.full(radii.size, np.nan)
    order = np.argsort(bins, kind="stable")
    bins_s, vals_s = bins[order], vals[order]
    starts = np.searchsorted(bins_s, np.arange(radii.size))
    ends = np.searchsorted(bins_s, np.arange(radii.size), side="right")
    for k in range(radii.size):
        v = vals_s[starts[k] : ends[k]]
        v = v[np.isfinite(v)]
        if v.size:
            prof[k] = np.max(v) if statistic == "max" else np.median(v)
    return radii, prof


def radial_profile_2d(
    im: np.ndarray,
    res: float = 1.0,
    rmin: float = 0.0,
    rmax: float | None = None,
    center: tuple[float, float] | None = None,
    *,
    med4q: bool = False,
    statistic: str = "median",
    nan_interp: bool = False,
) -> np.ndarray:
    """2-D radial-profile model (``profrad_tc`` ``p2d``).

    The 1-D profile is spline-interpolated at each pixel's radius; pixels
    outside ``[rmin, rmax]`` are 0 (as IDL) so subtracting the model leaves
    them untouched.  ``med4q`` instead returns the median of the four
    quadrant reflections of the image (``/med4q``).
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    if med4q:
        return _med4q(im)
    r = radial_grid(im.shape, center)
    if rmax is None:
        rmax = float(r.max())
    radii, prof = radial_profile(im, res, rmin, rmax, center, statistic=statistic)
    good = np.isfinite(prof)
    out = np.zeros_like(im)
    if good.sum() < 2:
        return out
    if nan_interp:
        prof_i = np.interp(radii, radii[good], prof[good])
        radii_i = radii
    else:
        prof_i, radii_i = prof[good], radii[good]
    spl = CubicSpline(radii_i, prof_i, bc_type="natural")
    inside = (r >= rmin) & (r <= rmax)
    vals = spl(np.clip(r[inside], radii_i[0], radii_i[-1]))
    out[inside] = vals
    if not nan_interp and (~good).any():
        # bins IDL marked NaN stay NaN in the model
        bad_r = radii[~good]
        rr = np.round(r / res) * res if res > 0 else r
        out[np.isin(rr, bad_r) & inside] = np.nan
    return out


def _med4q(im: np.ndarray) -> np.ndarray:
    ny, nx = im.shape
    xc, yc = image_center(im.shape)
    hx, hy = xc + 1, yc + 1
    q = np.full((4, hy, hx), np.nan)
    q[0] = im[:hy, :hx]
    q[1] = im[yc:, :hx][::-1][:hy]  # top-left reflected in y
    q[2] = im[yc:, xc:][::-1, ::-1][:hy, :hx]
    q[3] = im[:hy, xc:][:, ::-1][:, :hx]
    # IDL transposes quadrants 1 and 3 so that all four share the same orientation
    q[1] = q[1].T if q[1].shape[0] == q[1].shape[1] else q[1]
    q[3] = q[3].T if q[3].shape[0] == q[3].shape[1] else q[3]
    qm = np.nanmedian(q, axis=0)
    out = np.empty_like(im)
    out[:hy, :hx] = qm
    out[yc:, :hx] = (qm.T if qm.shape[0] == qm.shape[1] else qm)[::-1][: ny - yc]
    out[yc:, xc:] = qm[::-1, ::-1][: ny - yc, : nx - xc]
    out[:hy, xc:] = (qm.T if qm.shape[0] == qm.shape[1] else qm)[:, ::-1][:, : nx - xc]
    return out


def azimuthal_profile(
    im: np.ndarray,
    r1: float,
    r2: float,
    nbins: int = 36,
    center: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median vs. azimuth in the annulus ``[r1, r2]`` (``profang_tc``), plus a 2-D map.

    Returns ``(angles_deg, profile, map2d)``; the IDL version ignored ``r1``/``r2``
    (hard-coded 10–35) and left the last bin empty.
    """
    im = np.asarray(im, dtype=float)
    r = radial_grid(im.shape, center)
    from charis_dpp.geometry.grids import angle_grid

    th = angle_grid(im.shape, center, degrees=True) % 360.0
    edges = np.linspace(0, 360, nbins + 1)
    ang = 0.5 * (edges[:-1] + edges[1:])
    prof = np.full(nbins, np.nan)
    sel = (r >= r1) & (r <= r2)
    idx = np.clip(np.digitize(th, edges) - 1, 0, nbins - 1)
    for k in range(nbins):
        v = im[sel & (idx == k)]
        v = v[np.isfinite(v)]
        if v.size:
            prof[k] = np.median(v)
    good = np.isfinite(prof)
    map2d = np.zeros_like(im)
    if good.sum() >= 2:
        a_ext = np.concatenate([ang[good] - 360, ang[good], ang[good] + 360])
        p_ext = np.tile(prof[good], 3)
        map2d[sel] = np.interp(th[sel], a_ext, p_ext)
    return ang, prof, map2d


# --- bad-pixel repair --------------------------------------------------------------
def _box_sigma(box5: np.ndarray) -> float:
    srt = np.sort(box5.ravel())
    return float(np.min(srt[16:25] + srt[15:24] - srt[0:9] - srt[1:10]) / 2.0 / 1.540)


def _noise_rs(v: np.ndarray) -> float:
    s = np.sort(v)
    npix = s.size
    num = 2 * npix // 3 + 1
    noise = (
        np.min(
            s[num - 1 : npix]
            + s[num - 2 : npix - 1]
            - s[0 : npix - num + 1]
            - s[1 : npix - num + 2]
        )
        / 2.0
    )
    noise = noise / (1.9348 + ((num - 2) / npix - 0.66667) / 0.249891)
    if npix > 100:
        noise = noise / (1.0 - 0.6141 * np.exp(-0.5471 * np.log(npix)))
    else:
        noise = noise / (1.0 - 0.2223 * np.exp(-0.3294 * np.log(npix)))
    return float(noise)


def fixpix(im: np.ndarray, iterations: int = 3, return_mask: bool = False):
    """Iterative hot/cold pixel replacement on a 5×5 grid (``fixpix_rs.pro``).

    Returns the repaired image (and the mask of replaced pixels: +1 hot, -1 cold).
    Fixes the IDL sigma-buffer sizing bug (BUGS.md #267) and the inconsistent
    cold-pixel mask sign (#268).
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    imnew = im.copy()
    mask = np.zeros_like(im)
    medimg = np.nanmedian(im)
    sig_list = [
        _box_sigma(im[j - 2 : j + 3, i - 2 : i + 3])
        for i in range(2, nx - 2, 5)
        for j in range(2, ny - 2, 5)
        if np.isfinite(im[j - 2 : j + 3, i - 2 : i + 3]).all()
    ]
    if len(sig_list) < 4:
        return (imnew, mask) if return_mask else imnew
    sigma = np.array(sig_list)
    medsig = np.median(sigma)
    sigsig = _noise_rs(sigma)
    for it in range(1, iterations + 1):
        addon = ((it + 1) % 2) * 2
        for i in range(4 + addon, nx - 4, 5):
            for j in range(4 + addon, ny - 4, 5):
                box = imnew[j - 2 : j + 3, i - 2 : i + 3]
                if not np.isfinite(box).all():
                    continue
                for kind in ("hot", "cold"):
                    flat = np.argmax(box) if kind == "hot" else np.argmin(box)
                    py, px = divmod(int(flat), 5)
                    px += i - 2
                    py += j - 2
                    if px < 2 or py < 2 or px >= nx - 2 or py >= ny - 2:
                        continue
                    nb = imnew[py - 2 : py + 3, px - 2 : px + 3]
                    if not np.isfinite(nb).all():
                        continue
                    med8 = np.median(np.r_[nb[1, 1:4], nb[3, 1:4], nb[2, 1], nb[2, 3]])
                    med16 = np.median(np.r_[nb[0, :], nb[4, :], nb[1:4, 0], nb[1:4, 4]])
                    sig = _box_sigma(nb)
                    if sig > medsig + 2.0 * sigsig:
                        sig = max(medsig + 2.0 * sigsig, np.sqrt(5.0 + 0.210 * abs(med8)))
                    val = imnew[py, px]
                    if kind == "hot":
                        if (med8 + 2.0 * med16) / 3.0 - medimg > 2.0 * sig:
                            replace = val > 2.0 * med8 - medimg + 3.0 * sig
                        else:
                            replace = val - (med8 + 2.0 * med16) / 3.0 > 5.0 * sig
                        flag = 1
                    else:
                        if (med8 + 2.0 * med16) / 3.0 - medimg < -2.0 * sig:
                            replace = val < 2.0 * med8 - medimg - 3.0 * sig
                        else:
                            replace = val - (med8 + 2.0 * med16) / 3.0 < -5.0 * sig
                        flag = -1
                    if replace:
                        imnew[py, px] = med8
                        mask[py, px] = flag
    return (imnew, mask) if return_mask else imnew
