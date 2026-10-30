"""Image-quality helpers (``tools/imagequality/*.pro``)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.geometry.apertures import aperture_sum
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.filters import radial_profile, radial_profile_2d
from charis_dpp.io.fits import read_cube

__all__ = ["EncircledEnergy", "encircled_energy", "plane_fit_annulus"]


@dataclass
class EncircledEnergy:
    radii_px: np.ndarray
    ee: np.ndarray  # (nλ, nrad) fraction of the flux inside dpsf/2
    total_flux: np.ndarray
    profile_r: np.ndarray | None = None
    profile: np.ndarray | None = None


def encircled_energy(
    cube: np.ndarray | str | Path,
    *,
    position: tuple[float, float] | None = None,
    dpsf: int = 41,
    channel: int | None = None,
    prad: bool = False,
    background: float = 0.0,
    plot: str | Path | None = None,
) -> EncircledEnergy:
    """Encircled-energy curves per channel (``charis_calc_ee.pro``, with the missing
    ``profradxy`` replaced by the radial profile about the source, #259)."""
    data = read_cube(cube).data if isinstance(cube, str | Path) else np.asarray(cube, dtype=float)
    if data.ndim == 2:
        data = data[None]
    nlam = data.shape[0]
    dpsf = min(dpsf, data.shape[-1])
    x0, y0 = image_center(data.shape) if position is None else position
    radii = np.arange(dpsf) / 2.0
    ee = np.zeros((nlam, dpsf))
    tot = np.zeros(nlam)
    for il in range(nlam):
        im = data[il] - background
        tot[il] = aperture_sum(im, x0, y0, dpsf / 2.0)
        for j, r in enumerate(radii):
            ee[il, j] = aperture_sum(im, x0, y0, r) / tot[il] if r > 0 else 0.0
    il = nlam // 2 if channel is None else channel
    im = data[il] - background
    if prad:
        im = im - radial_profile_2d(im, 1.0, 1.0, None, None)
    rr, pr = radial_profile(im, 1.0, 1.0, dpsf / 2.0, (x0, y0))
    core = radial_grid(im.shape, (x0, y0)) <= dpsf / 4.0
    pr = pr / np.nanmax(im[core])
    res = EncircledEnergy(radii, ee, tot, rr, pr)
    if plot is not None:
        from charis_dpp.plotting import finish, plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(radii, ee[il], "-", label=f"encircled energy, channel {il}")
        ax.plot(rr, pr, ":", label="normalised radial profile")
        ax.set_xlabel("radius [px]")
        ax.legend(frameon=False)
        finish(fig, Path(plot))
    return res


def plane_fit_annulus(
    image: np.ndarray, source: tuple[float, float], rad: float, *, center=None
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a plane to the annular background around a source and return
    ``(background-subtracted image, fitted plane)`` (``charis_planefit_ex.pro``, generalised)."""
    im = np.asarray(image, dtype=float)
    center = image_center(im.shape) if center is None else center
    rs = radial_grid(im.shape, source)
    rc = radial_grid(im.shape, center)
    sep = float(np.hypot(source[0] - center[0], source[1] - center[1]))
    src = rs <= 2 * rad
    bkg = (rs <= 4 * rad) & ~src & (np.abs(rc - sep) < 3 * rad) & np.isfinite(im)
    yy, xx = np.nonzero(bkg)
    a = np.column_stack([np.ones(xx.size), xx, yy])
    coeff, *_ = np.linalg.lstsq(a, im[bkg], rcond=None)
    gy, gx = np.mgrid[0 : im.shape[0], 0 : im.shape[1]]
    plane = coeff[0] + coeff[1] * gx + coeff[2] * gy
    plane = np.where(src | bkg, plane, 0.0)
    return im - plane, plane
