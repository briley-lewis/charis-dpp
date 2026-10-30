"""Simple classical SDI (``charis_ssdi.pro``).

For each cube and channel ``il``, all other channels are speckle-aligned to
``il`` and their median (optionally least-squares scaled) is subtracted.
Fixes vs. IDL: works on any number of channels (#156), the reference median
excludes the channel itself (#157), the north-up result is actually written (#155).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import align_speckles
from charis_dpp.geometry.grids import image_center
from charis_dpp.image.filters import idl_filter_image
from charis_dpp.image.stats import combine_frames
from charis_dpp.image.wcs import northup
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["ssdi"]

SOURCES = {"aloci": ("proc", "_alocisub"), "klip": ("proc", "_adisub"), "rsub": ("rsub", "rsub")}


def ssdi(
    cfg: DatasetConfig | str | Path,
    *,
    source: str = "aloci",
    weight: bool = False,
    filt: bool = False,
    northup_final: bool = False,
    prefix: str = "n",
    outfile: str = "ssdi.fits",
) -> Path:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    sub, suffix = SOURCES[source]
    nums = cfg.frames("science")
    wav = res.wavelengths_um
    center = res.star_center or None
    frames, h = [], None
    for n in nums:
        c = read_cube(getattr(dirs, sub) / frame_filename(n, suffix, prefix))
        cube = c.data
        nlam = cube.shape[0]
        out = np.empty_like(cube)
        for il in range(nlam):
            aligned = align_speckles(cube, wav, il, center or image_center(c.shape))
            others = np.delete(aligned, il, axis=0)
            ref = np.nanmedian(others, axis=0)
            w = 1.0
            if weight:
                w = np.nansum(cube[il] * ref) / np.nansum(ref * ref)
            d = cube[il] - w * ref
            if filt:
                d = d - idl_filter_image(d, median=15)
            out[il] = d
        c.h1["SSDI"] = (1, "simple SDI applied")
        if source == "rsub" and northup_final:
            out, _ = northup(out, c.h1, res.instrument.north_offset_deg, center=center)
        write_cube(getattr(dirs, sub) / frame_filename(n, "_ssdisub", prefix), out, c.h0, c.h1)
        frames.append(out)
        h = h or (c.h0, c.h1)
    combined = combine_frames(np.stack(frames), "median", axis=0)
    outp = dirs.proc / outfile
    write_cube(outp, combined, *h)
    write_cube(outp.with_name(outp.stem + "_collapsed.fits"), np.nanmedian(combined, axis=0), *h)
    log.info("ssdi: wrote %s", outp)
    return outp
