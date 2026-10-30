"""PSF quality diagnostic (``charis_check_psfquality.pro`` / ``register_cube,/checkquality``).

For each registered cube the median satellite-spot flux (aperture of radius
``fwhm/2``) is divided by the median halo brightness at 0.2–0.4″; the ratio
per frame is written to ``reduc/reg/PSFquality.txt`` and plotted.  All
channels are used when ``all_slices`` (the IDL version looked at channel 0
only, BUGS.md #47).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.image.filters import radial_profile_2d
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import get_psf_center, get_sat_positions, read_cube
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["check_psf_quality"]


def check_psf_quality(
    cfg: DatasetConfig | str | Path,
    *,
    suffix: str = "reg",
    prefix: str = "n",
    all_slices: bool = False,
    halo_arcsec: tuple[float, float] = (0.2, 0.4),
    bin_size: float = 0.5,
    plot: bool = True,
) -> np.ndarray:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    nums = cfg.frames("science")
    pixscale = res.instrument.pixscale_arcsec
    ratios = np.full(len(nums), np.nan)
    for i, n in enumerate(nums):
        p = dirs.reg / frame_filename(n, suffix, prefix)
        if not p.exists():
            continue
        c = read_cube(p)
        sats = get_sat_positions(c.h1, c.nlam)
        xc, yc = get_psf_center(c.h1) or (c.shape[1] // 2, c.shape[0] // 2)
        chans = range(c.nlam) if all_slices else [0]
        vals = []
        for il in chans:
            im = c.data[il]
            ap = 0.5 * res.fwhm_px[il]
            fl = [
                aper_idl(im, sx, sy, ap, setskyval=0.0)[0] for sx, sy in sats[il] if np.isfinite(sx)
            ]
            if not fl:
                continue
            halo = radial_profile_2d(im, 1.0, 1.0, min(im.shape) // 2 - 1, (xc, yc))
            r = radial_grid(im.shape, (xc, yc))
            sel = (r >= halo_arcsec[0] / pixscale) & (r <= halo_arcsec[1] / pixscale)
            hv = np.nanmedian(halo[sel]) if sel.any() else np.nan
            if np.isfinite(hv) and hv != 0:
                vals.append(np.median(fl) / hv)
        if vals:
            ratios[i] = float(np.median(vals))
    out = dirs.reg / "PSFquality.txt"
    out.write_text("".join(f"{n:5d} {r:12.5f}\n" for n, r in zip(nums, ratios, strict=True)))
    if plot:
        from charis_dpp.plotting import finish, plt

        fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
        ax[0].plot(nums, ratios, "o")
        ax[0].set_xlabel("frame")
        ax[0].set_ylabel("sat / halo")
        good = np.isfinite(ratios)
        if good.any():
            ax[1].hist(ratios[good], bins=max(int(np.ptp(ratios[good]) / bin_size), 1))
        ax[1].set_xlabel("sat / halo")
        finish(fig, dirs.reg / "psfquality.png")
    log.info("PSF quality: median sat/halo %.3f (%s)", np.nanmedian(ratios), out)
    return ratios
