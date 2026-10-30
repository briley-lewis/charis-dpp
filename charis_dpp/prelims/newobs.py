"""``charis newobs`` — create the dataset directory layout and configuration (``charis_newobs.pro``)."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from charis_dpp.config.dataset import DatasetConfig, canonical_band
from charis_dpp.config.instrument import to_mjd
from charis_dpp.io.catalogs import StarInfo, query_star
from charis_dpp.prelims.setupdir import setup_dirs

log = logging.getLogger(__name__)

__all__ = ["newobs"]


def newobs(
    date: str,
    object_name: str,
    band: str,
    *,
    workdir: str | Path = ".",
    star: StarInfo | None = None,
    astrogrid_lambda_d: float | None = None,
    modulation_nm: float | None = None,
    move_raw: bool = True,
) -> Path:
    """Set up ``workdir`` for a new observation and write ``<Object>_<band>.toml`` (+ legacy ``.info``).

    ``star`` bypasses the SIMBAD/VizieR lookup (offline use, tests).  Raw
    ``CRSA*cube.fits`` files in ``workdir`` are moved into ``data/raw/``.
    Returns the config path.
    """
    workdir = Path(workdir)
    dirs = setup_dirs(workdir)
    if move_raw:
        for f in sorted(workdir.glob("CRSA*cube.fits")):
            shutil.move(str(f), dirs.raw / f.name)
    b = canonical_band(band)
    if b is None:
        raise ValueError(f"unknown band {band!r} (use lowres/low/J/H/K)")
    if star is None:
        log.info("resolving %s with SIMBAD/VizieR", object_name)
        star = query_star(object_name)
    mjd = to_mjd(date)
    cfg = DatasetConfig.new(
        target=star.to_config(),
        observation={
            "date": date if "-" in date else f"{date[:4]}-{date[4:6]}-{date[6:]}",
            "mjd": mjd,
            "band": b,
            "astrogrid_lambda_d": astrogrid_lambda_d,
            "modulation_nm": modulation_nm,
            "frames": {"science": "", "sky": "", "ladder": ""},
        },
    )
    stem = f"{object_name.replace(' ', '')}_{'low' if b == 'lowres' else b}"
    path = cfg.save(workdir / f"{stem}.toml")
    log.info("wrote %s and %s", path, path.with_suffix(".info"))
    return path
