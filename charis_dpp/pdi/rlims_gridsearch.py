"""Grid search over the unresolved-polarisation annulus radii
(``charis_pdi_rlims_gridsearch.pro``; the dataset is a parameter, #227)."""

from __future__ import annotations

import logging
from pathlib import Path

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.pdi.pdi import pdi
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["rlims_gridsearch"]


def rlims_gridsearch(
    cfg: DatasetConfig | str | Path,
    rmin_vals=(0, 2, 4, 6, 8, 10),
    rmax_vals=(5, 10, 15, 20, 25, 30),
    *,
    logfile: str = "pdi_rlims_gridsearch.log",
    **pdi_kwargs,
) -> list[tuple[int, int, str]]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    dirs = Dirs(cfg.workdir)
    runs = []
    lines = []
    for rmin in rmin_vals:
        for rmax in rmax_vals:
            if rmin >= rmax:
                continue
            out = f"pdi_rmin{int(rmin):02d}_rmax{int(rmax):02d}"
            kw = dict(meanadd=True, bg_sub=True, unres_pol_sub=True, plot=False)
            kw.update(pdi_kwargs)
            pdi(cfg, unres_pol_rlims=(rmin, rmax), outfile=out, **kw)
            runs.append((int(rmin), int(rmax), out))
            lines.append(f"{rmin} {rmax} {out}")
            log.info("grid search: rmin=%s rmax=%s -> %s", rmin, rmax, out)
    (dirs.root / logfile).write_text("\n".join(lines) + "\n")
    return runs
