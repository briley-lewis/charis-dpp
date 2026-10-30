"""Shared machinery for the PSF-subtraction drivers.

The four IDL codes (``charis_adialoci``, ``charis_sdialoci``, ``charis_adiklip``,
``charis_rdiklip``) duplicated ~500 lines of setup, ring-caching and output
code each; it lives here once.  The IDL scratch-file protocol
(``indices_*.dat`` / ``values_*.dat``) is replaced by in-memory per-channel
stacks and boolean masks.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.io.fits import Header

from charis_dpp.astro.time_coords import rotation_rate
from charis_dpp.config.dataset import DatasetConfig, Resolved
from charis_dpp.geometry.grids import angle_grid, image_center, radial_grid
from charis_dpp.image.stats import combine_frames
from charis_dpp.image.wcs import northup
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = [
    "PsfSubSetup",
    "prepare",
    "load_channel",
    "ring_edges",
    "loci_geometry",
    "write_results",
    "HEADER_KEYS",
]

HEADER_KEYS = {
    "aloci": [
        "loci_nfw",
        "loci_na",
        "loci_geo",
        "loci_drs",
        "loci_rmi",
        "loci_rma",
        "loci_nre",
        "svd",
        "pixmaskv",
        "zero",
        "meanadd",
        "rsub",
    ],
    "klip": [
        "klip_pca",
        "klip_nfw",
        "klip_drs",
        "klip_rmi",
        "klip_rma",
        "klip_nzo",
        "meansub",
        "zero",
        "meanadd",
        "rsub",
    ],
}


@dataclass
class PsfSubSetup:
    cfg: DatasetConfig
    res: Resolved
    dirs: Dirs
    algorithm: str
    nums: list[int]
    files: list[Path]
    files_fwd: list[Path] | None
    out_files: list[Path]
    outfile: Path
    pa: np.ndarray
    ha: np.ndarray
    dtpose: np.ndarray  # smear during one exposure, degrees
    rsubval: int
    center: tuple[float, float]
    shape: tuple[int, int]
    params: dict = field(default_factory=dict)
    headers: list[tuple[Header, Header]] = field(default_factory=list)
    norot: bool = False
    meanadd: bool = False
    zero: bool = False

    @property
    def nfiles(self) -> int:
        return len(self.files)

    @property
    def nlam(self) -> int:
        return self.res.nlam

    @property
    def wav(self) -> np.ndarray:
        return self.res.wavelengths_um

    @property
    def fwhm(self) -> np.ndarray:
        return self.res.fwhm_px


def prepare(
    cfg: DatasetConfig | str | Path,
    *,
    algorithm: str,
    rsubval: int | None = None,
    fc: bool = False,
    fwdmod: bool = False,
    postsdi: bool = False,
    postadi: bool = False,
    prefix: str = "n",
    suffix: str | None = None,
    outfile: str | None = None,
    out_suffix: str = "_alocisub",
    norot: bool = False,
    meanadd: bool | None = None,
    zero: bool | None = None,
) -> PsfSubSetup:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    red = cfg.reduction
    rsubval = int(rsubval if rsubval is not None else bool(red.get("rsub", False)))
    nums = cfg.frames("science")
    if not nums:
        raise ValueError("no science frames in observation.frames")
    if postsdi or postadi:
        datadir = dirs.proc
        suffix = suffix or ("_sdialocisub" if postsdi else "_alocisub")
    elif rsubval > 0:
        datadir = dirs.rsub
        suffix = suffix or "rsub"
    else:
        datadir = dirs.reg
        if suffix is None:
            suffix = (
                "reg_cal"
                if (dirs.reg / frame_filename(nums[0], "reg_cal", prefix)).exists()
                else "reg"
            )
    if fc and not fwdmod:
        files = [datadir / frame_filename(n, suffix + "_fc", prefix) for n in nums]
        files_fwd = None
    else:
        files = [datadir / frame_filename(n, suffix, prefix) for n in nums]
        files_fwd = (
            [datadir / frame_filename(n, suffix + "_fc", prefix) for n in nums] if fwdmod else None
        )
    for p in files + (files_fwd or []):
        if not p.exists():
            raise FileNotFoundError(f"{p} not found")
    tag = "_fc" if fc else ""
    out_files = [dirs.proc / frame_filename(n, out_suffix + tag, prefix) for n in nums]
    outfile = dirs.proc / (outfile or f"final{tag}.fits")
    dirs.proc.mkdir(parents=True, exist_ok=True)

    log_ = ReducLog.load(dirs.reduclog)
    log_.require(nums)
    pa = log_.pa(nums)
    ha = log_.ha(nums)
    first = read_cube(files[0])
    if "TOT_ROT" not in first.h0 and "CD1_1" not in first.h1 and "CD1_1" not in first.h0:
        raise ValueError(
            f"{files[0]}: no WCS rotation information (TOT_ROT / CD matrix) in the headers"
        )
    inst = res.instrument
    dec = float(cfg.target.get("dec_deg", first.get("DEC_DEG", 0.0)))
    tint = np.array(
        [
            float(read_cube(p).get("EXP1TIME", inst.read_time_s))
            * float(read_cube(p).get("COADDS", 1))
            for p in files
        ]
    )
    dtpose = np.abs(rotation_rate(ha, dec, inst.latitude_deg)) * tint
    center = res.star_center or image_center(first.shape)
    setup = PsfSubSetup(
        cfg=cfg,
        res=res,
        dirs=dirs,
        algorithm=algorithm,
        nums=nums,
        files=files,
        files_fwd=files_fwd,
        out_files=out_files,
        outfile=outfile,
        pa=pa,
        ha=ha,
        dtpose=dtpose,
        rsubval=rsubval,
        center=(float(center[0]), float(center[1])),
        shape=first.shape,
        norot=norot,
        meanadd=bool(red.get("meanadd", False) if meanadd is None else meanadd),
        zero=bool(red.get("zero", False) if zero is None else zero),
    )
    setup.headers = [(read_cube(p).h0, read_cube(p).h1) for p in files]
    return setup


def load_channel(setup: PsfSubSetup, il: int, *, model: bool = False) -> np.ndarray:
    """``(nf, ny, nx)`` stack of channel ``il`` from the data (or forward-model) cubes."""
    paths = setup.files_fwd if model else setup.files
    out = np.empty((len(paths),) + setup.shape, dtype=float)
    for i, p in enumerate(paths):
        with fits.open(p, memmap=True) as hdul:
            hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None else hdul[0]
            out[i] = np.asarray(hdu.data[il], dtype=float)
    return out


def ring_edges(rmin: float, rmax: float, drsub) -> list[tuple[float, float]]:
    """Subtraction annuli ``[(ri, rf), ...]``.

    ``drsub`` is a scalar width, or the IDL 4-vector ``[w_in, w_out, r_break, r_scale]``
    giving ``dr(r) = w_in + (w_out - w_in)·(0.5 + atan((r - r_break)/r_scale)/π)``.
    (The IDL ``n_elements(drsub0 eq 1)`` test made the vector form unreachable, BUGS.md #127.)
    """
    d = np.atleast_1d(np.asarray(drsub, dtype=float))
    edges = []
    if d.size == 1:
        w = float(d[0])
        n = int(np.ceil((rmax - rmin) / w))
        for i in range(max(n, 1)):
            ri = rmin + i * w
            edges.append((ri, min(ri + w, rmax)))
        return edges
    if d.size != 4:
        raise ValueError("drsub must be a scalar or [w_in, w_out, r_break, r_scale]")
    r = float(rmin)
    while r < rmax:
        dr = (0.5 + np.arctan((r - d[2]) / d[3]) / np.pi) * (d[1] - d[0]) + d[0]
        edges.append((r, min(r + dr, rmax)))
        r += dr
    return edges


def loci_geometry(
    ri: float, dr: float, fwhm: float, na, geom: float, rim_max: float
) -> tuple[int, float, float]:
    """Lafrenière optimisation-zone geometry for a ring starting at ``ri``: ``(nt, r1opt, r2opt)``.

    ``na`` is the optimisation area in PSF cores (scalar or the IDL 4-vector
    arctan profile in radius).
    """
    n = np.atleast_1d(np.asarray(na, dtype=float))
    r = ri + dr / 2.0
    if n.size == 1:
        area = float(n[0]) * np.pi * (fwhm / 2.0) ** 2
    else:
        area = (
            ((0.5 + np.arctan((r - n[2]) / n[3]) / np.pi) * (n[1] - n[0]) + n[0])
            * np.pi
            * (fwhm / 2.0) ** 2
        )
    dropt = np.sqrt(geom * area)
    r1opt = ri
    nt = max(int(round((2 * np.pi * (r1opt + dropt / 2.0) * dropt) / area)), 1)
    dropt = np.sqrt(r1opt**2 + (nt * area) / np.pi) - r1opt
    r2opt = r1opt + dropt
    if r2opt > rim_max:
        r2opt = rim_max
        dropt = r2opt - r1opt
        nt = max(int(round((2 * np.pi * (r1opt + dropt / 2.0) * dropt) / area)), 1)
    return nt, r1opt, r2opt


def zone_masks(shape, center, ri, rf, r1opt, r2opt, nt, it, *, pixmask: bool):
    """Boolean ``(sub, opt)`` masks for sector ``it`` of ``nt`` (IDL sector angles from +x, counter-clockwise)."""
    r = radial_grid(shape, center)
    ang = angle_grid(shape, center) % (2 * np.pi)
    dt = 2 * np.pi / nt
    sector = (ang >= it * dt) & (ang < (it + 1) * dt)
    opt = (r >= r1opt) & (r < r2opt) & sector
    sub = (r >= ri) & (r < rf) & sector
    if pixmask:
        opt = opt & (r >= r1opt + (rf - ri))
    return sub, opt


def write_results(
    setup: PsfSubSetup, resid: np.ndarray, *, params: dict, history: str
) -> tuple[Path, Path]:
    """Write per-frame residual cubes (north-up unless ``norot``), the combined cube and its collapse.

    ``resid`` is ``(nf, nλ, ny, nx)``.  Header parameter cards are written with
    the IDL names (truncated to 8 characters) so ``partialsub`` and the
    detection/contrast tools can read them back.
    """
    inst = setup.res.instrument
    frames = []
    h1_first = None
    for i, (h0, h1) in enumerate(setup.headers):
        h0, h1 = h0.copy(), h1.copy()
        data = resid[i]
        if not setup.norot:
            data, rot = northup(data, h1, inst.north_offset_deg, center=setup.center, frame_index=i)
            h1.add_history(f"A rotation of {rot:.3f} degrees was applied to the image.")
        else:
            h1["ROTANG"] = (0.0, "no north-up rotation applied (norot)")
        for k, v in params.items():
            h1[k[:8].upper()] = v
        h1["RSUB"] = setup.rsubval
        h1["PSFALGO"] = (setup.algorithm, "PSF subtraction algorithm")
        h0.add_history(history)
        write_cube(setup.out_files[i], data, h0, h1)
        frames.append(data)
        if h1_first is None:
            h0_first, h1_first = h0, h1
    stack = np.stack(frames)
    combined = combine_frames(stack, "mean" if setup.meanadd else "median", axis=0)
    write_cube(setup.outfile, combined, h0_first, h1_first)
    collapsed = combine_frames(combined, "mean" if setup.meanadd else "median", axis=0)
    col = setup.outfile.with_name(setup.outfile.stem + "_collapsed.fits")
    write_cube(col, collapsed, h0_first, h1_first)
    log.info("wrote %s and %s", setup.outfile, col)
    return setup.outfile, col


def outlier_mask(optreg: np.ndarray, clip: float = 7.0) -> np.ndarray:
    """Pixels (rows of an ``(npix, nf)`` matrix) finite in every frame and within ``clip``× the robust scale."""
    z = optreg.copy()
    with np.errstate(invalid="ignore", divide="ignore"):
        col = np.nanmedian(np.abs(z), axis=0)
        z = z / np.where(col > 0, col, 1.0)[None, :]
        row = np.nanmedian(np.abs(z), axis=1)
        z = z / np.where(row > 0, row, 1.0)[:, None]
    ok = np.isfinite(z) & (np.abs(z) < clip)
    return ok.all(axis=1)
