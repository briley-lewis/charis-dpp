"""Header preparation (``charis_imprep.pro``).

For every raw cube in ``data/raw``: normalise the header (HA, LST, exposure
time, coadds, site, PA, wavelength axis), classify it as science/sky/bad,
optionally repair bad pixels, apply the post-NBS parity correction when the
instrument epoch calls for it, and write ``reduc/prep/n####e.fits``.  Writes
``obslog.txt`` and records the science/sky frame lists in the dataset config.

Fixes vs. IDL: RA/DEC are always read (the IDL version crashed when the header
already had LST, BUGS.md #11); HA is computed consistently (#12/#13); only
``CRSA*_cube.fits`` files are numbered, so stray FITS files cannot shift the
frame numbers (#14); the obslog header matches its columns (#15).
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits

from charis_dpp.astro.time_coords import hour_angle, lst_hours, precess, ten
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.config.instrument import Instrument, instrument_for
from charis_dpp.image.filters import fixpix as _fixpix
from charis_dpp.image.wcs import nbs_parity_correct
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import header_get, write_cube
from charis_dpp.io.obslog import ObsEntry, ObsLog
from charis_dpp.prelims.health import classify_frame
from charis_dpp.prelims.setupdir import Dirs, setup_dirs

log = logging.getLogger(__name__)

__all__ = ["imprep", "RAW_EXT3_KEYS"]

RAW_EXT3_KEYS = ("D_IMRANG", "RET-ANG1", "MJD-STR", "MJD-END")


def _load_cfg(cfg: DatasetConfig | str | Path) -> DatasetConfig:
    return cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)


def imprep(
    cfg: DatasetConfig | str | Path,
    *,
    fixpix: bool = False,
    pdi: bool = False,
    manual_exptime: float | None = None,
    manual_coadds: int | None = None,
    prefix: str = "n",
    raw_pattern: str = "CRSA*_cube.fits",
) -> ObsLog:
    cfg = _load_cfg(cfg)
    dirs: Dirs = setup_dirs(cfg.workdir)
    for f in sorted(cfg.workdir.glob("CRSA*fits")):
        shutil.move(str(f), dirs.raw / f.name)
    files = sorted(dirs.raw.glob(raw_pattern)) or sorted(dirs.raw.glob("*.fit*"))
    if not files:
        raise FileNotFoundError(f"no raw cubes matching {raw_pattern} in {dirs.raw}")
    log.info("imprep: %d raw cubes in %s", len(files), dirs.raw)

    obslog = ObsLog(path=dirs.obslog)
    inst_cache: dict[str, Instrument] = {}
    for i, f in enumerate(files):
        with fits.open(f, memmap=False) as hdul:
            h0 = hdul[0].header.copy()
            if len(hdul) > 1 and hdul[1].data is not None:
                cube = np.asarray(hdul[1].data, dtype=float)
                h1 = hdul[1].header.copy()
            else:
                cube = np.asarray(hdul[0].data, dtype=float)
                h1 = h0.copy()
            if len(hdul) > 3:
                h3 = hdul[3].header
                for k in RAW_EXT3_KEYS:
                    if k in h3 and k not in h0:
                        h0[k] = (h3[k], h3.comments[k])
        mjd = h0.get("MJD", h1.get("MJD"))
        if mjd is None:
            raise KeyError(f"{f.name}: MJD missing from header; cannot prepare the cube")
        inst = inst_cache.get(f"{mjd:.3f}") or instrument_for(
            float(mjd), cfg.data.get("instrument") or None
        )
        inst_cache[f"{mjd:.3f}"] = inst

        exp1time = manual_exptime if manual_exptime is not None else inst.read_time_s
        if manual_coadds is not None:
            coadds = int(manual_coadds)
        else:
            first = header_get(h1, "FIRSTRD", header_get(h0, "FIRSTRD", 1))
            last = header_get(h1, "LASTRD", header_get(h0, "LASTRD", 1))
            coadds = int(last - first + 1)

        # --- time / pointing ---------------------------------------------------
        ra = 15.0 * ten(header_get(h0, "RA"))
        dec = ten(header_get(h0, "DEC"))
        lst = None
        for key in ("ST", "LST"):
            if key in h0:
                lst = ten(h0[key])
        if lst is None:
            lst = float(lst_hours(float(mjd), inst.longitude_deg))
        try:
            from astropy.time import Time

            epoch = float(Time(float(mjd), format="mjd").jyear)
        except Exception:  # pragma: no cover
            epoch = 2000.0
        rap, decp = precess(ra, dec, 2000.0, epoch)
        ha = float(hour_angle(lst, rap))
        h0["HA"] = (ha, "hour angle [hours]")
        h0["LST"] = (lst, "local sidereal time [hours]")
        h0["EXP1TIME"] = (exp1time, "single read time [s]")
        h0["COADDS"] = (coadds, "number of reads")
        h0["TELESCOP"] = "Subaru"
        h0["LAT"] = (inst.latitude_deg, "site latitude [deg]")
        h0["LNG"] = (inst.longitude_deg, "site longitude [deg]")
        pa0 = float(header_get(h0, "PARANG", header_get(h1, "PARANG", 0.0)))
        h0["PA"] = (pa0, "parallactic angle [deg]")
        h0["RA_DEG"] = (ra, "RA J2000 [deg]")
        h0["DEC_DEG"] = (dec, "Dec J2000 [deg]")
        for k in ("HA", "LST", "EXP1TIME", "COADDS", "PA", "RA_DEG", "DEC_DEG"):
            h1[k] = (h0[k], h0.comments[k])
        h1["CRPIX3"] = 1
        h1["CUNIT3"] = "microns"
        for k, v, c in inst.to_header_cards():
            h0[k] = (v, c)
        h0.add_history(f"charis imprep: {f.name} -> {frame_filename(i + 1, 'e', prefix)}")

        health = classify_frame(cube, inst, pdi=pdi)
        if fixpix:
            cube = np.stack([_fixpix(sl, 10) for sl in cube])
        if inst.post_nbs:
            cube = nbs_parity_correct(cube, h0, h1, float(inst.get("nbs.rotzp_deg")), parang=pa0)
        else:
            h0["PARITY"] = (0, "0: as delivered (pre-NBS); 1: x-flipped by the NBS correction")

        out = frame_filename(i + 1, "e", prefix)
        write_cube(dirs.prep / out, cube, h0, h1)
        obslog.entries.append(
            ObsEntry(
                f.name,
                out,
                str(h0.get("FILTNAME", h0.get("CAL_BAND", "?"))).replace(" ", "_"),
                exp1time * coadds,
                pa0,
                health,
            )
        )
        log.info(
            "%s -> %s  exptime %.1fs  PA %.2f  %s", f.name, out, exp1time * coadds, pa0, health
        )
    obslog.save()
    cfg.set_frames("science", obslog.frames("science"))
    cfg.set_frames("sky", obslog.frames("sky"))
    if cfg.geometry.get("frame_shape") is None:
        cfg.set("geometry.frame_shape", [int(cube.shape[-2]), int(cube.shape[-1])])
    if cfg.observation.get("mjd") is None:
        cfg.set("observation.mjd", float(mjd))
    if cfg.path is not None:
        cfg.save()
    return obslog
