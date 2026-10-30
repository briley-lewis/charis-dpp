"""Astrometry from a grid of starting positions (``charis_est_astrom_grid.pro``):
extract at every (x, y) in a box, keep the centroid results, and report the
mean and SNR-weighted position.  The IDL defaults were one target's values
(#119); here the box is required.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.extract_spectrum import extract_spectrum

log = logging.getLogger(__name__)

__all__ = ["astrometry_grid"]


def astrometry_grid(
    cfg: DatasetConfig | str | Path,
    datacube: str | Path,
    xrange: tuple[float, float],
    yrange: tuple[float, float],
    *,
    deltapos: float = 1.0,
    fitrange: list[int] | None = None,
    bklim: float = 3.0,
    bklir: float = 1.0,
    throughput: str | Path | None = None,
) -> dict:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    dirs = Dirs(cfg.workdir)
    rows = []
    for x0 in np.arange(xrange[0], xrange[1] + 1e-9, deltapos):
        for y0 in np.arange(yrange[0], yrange[1] + 1e-9, deltapos):
            try:
                r = extract_spectrum(
                    cfg,
                    datacube,
                    (float(x0), float(y0)),
                    fitbckgd=True,
                    fitrange=fitrange,
                    bklim=bklim,
                    bklir=bklir,
                    throughput=throughput,
                    plot=False,
                    outputspec="astrom_tmp.dat",
                )
            except RuntimeError:
                continue
            rows.append((x0, y0, r.position[0], r.position[1], r.snr_collapsed))
    if not rows:
        raise RuntimeError("no successful extractions in the grid")
    arr = np.array(rows)
    np.savetxt(dirs.proc / "astrom_output.txt", arr, fmt="%.4f", header="x0 y0 xfit yfit snr")
    xf, yf, snr = arr[:, 2], arr[:, 3], arr[:, 4]
    w = np.clip(snr, 0, None) ** 2
    out = {
        "x_mean": float(np.mean(xf)),
        "y_mean": float(np.mean(yf)),
        "x_std": float(np.std(xf)),
        "y_std": float(np.std(yf)),
        "x_weighted": float(np.sum(w * xf) / np.sum(w)) if w.sum() > 0 else np.nan,
        "y_weighted": float(np.sum(w * yf) / np.sum(w)) if w.sum() > 0 else np.nan,
        "n": int(arr.shape[0]),
    }
    (dirs.proc / "avgastrom.txt").write_text(
        f"{out['x_mean']:.4f} {out['y_mean']:.4f} {out['x_std']:.4f} {out['y_std']:.4f}\n"
    )
    (dirs.proc / "weighted_astrom.txt").write_text(
        f"{out['x_weighted']:.4f} {out['y_weighted']:.4f}\n"
    )
    (dirs.proc / "astrom_tmp.dat").unlink(missing_ok=True)
    log.info(
        "astrometry grid: (%.3f ± %.3f, %.3f ± %.3f) from %d fits",
        out["x_mean"],
        out["x_std"],
        out["y_mean"],
        out["y_std"],
        out["n"],
    )
    return out
