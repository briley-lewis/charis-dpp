"""Double-difference polarimetry (``charis_pdi.pro``, fork version).

For each complete HWP cycle (positions 0/22.5/45/67.5°): single sums and
differences of the two beams give ``I_Q, I_U, Q, U``; the Mueller-matrix
model of the train gives, per channel, the double-difference row vector
``(m_I, m_Q, m_U)`` for each cycle, from which the instrumental polarisation
(``Q −= m_I·I_Q``) and the Q/U crosstalk (2×2 solve) are removed.  Cycles are
derotated north-up, optionally background- and unresolved-stellar-
polarisation-subtracted, and combined.

Outputs under ``reduc/proc/``: ``<out>_q/u/pi/aolp/qphi/uphi.fits`` cubes and
their ``_collapsed`` images, ``coverage_map.fits`` and optional per-cycle
products.

Fixes vs. IDL: no ``where() → -1`` pixel corruption (#206); an incomplete
cycle is skipped, not silently used (#207/#208); frames are matched by
number, not log row (#209); header cards are valid (#210/#211); telluric
channels and the φ grid come from the data (#212/#213); temporary files are
removed with the standard library (#214).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.stats import combine_frames, robust_poly_fit
from charis_dpp.image.wcs import northup
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import get_sat_positions, header_get, read_cube, write_cube
from charis_dpp.pdi.hwp_match import HwpTable
from charis_dpp.pdi.mueller.matrix import MuellerTerms, mueller_rows_per_channel
from charis_dpp.pdi.stokes import (
    adjust_pol_axis,
    average_angle,
    azimuthal_stokes,
    polarisation_angle_grid,
)
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["pdi", "PdiResult"]


@dataclass
class PdiResult:
    q: np.ndarray
    u: np.ndarray
    qphi: np.ndarray
    uphi: np.ndarray
    outdir: Path
    ncycles: int
    files: dict[str, Path] = field(default_factory=dict)
    unres_q: np.ndarray | None = None  # (nλ, ncycles)
    unres_u: np.ndarray | None = None


def _collapse(cube: np.ndarray, how: str) -> np.ndarray:
    if how == "median":
        return np.nanmedian(cube, axis=0)
    if how == "mean":
        return np.nanmean(cube, axis=0)
    if how == "rmean":
        return combine_frames(cube, "mean", axis=0)
    raise ValueError("collapse_type must be median, mean or rmean")


def _write(path: Path, data, h0: Header, h1: Header) -> Path:
    write_cube(path, data, h0, h1)
    return path


def pdi(
    cfg: DatasetConfig | str | Path,
    *,
    hwp_info: str | Path | HwpTable = "hwp_info.txt",
    outfile: str = "pdi",
    prefix: str = "n",
    left_suffix: str | None = None,
    right_suffix: str | None = None,
    meanadd: bool = False,
    collapse_type: str | None = None,
    angoffset: float | None = None,
    maskspots: bool = False,
    mm_ip_corr: bool = True,
    mm_crosstalk_corr: bool = True,
    first_order_ip_corr: bool = False,
    foip_rlims: tuple[float, float] | None = None,
    bg_sub: bool = False,
    bg_rlims: tuple[float, float] | None = None,
    unres_pol_sub: bool = False,
    unres_pol_rlims: tuple[float, float] | None = None,
    avg_unres_pol_coeffs: bool = False,
    collapse_cycles: bool = False,
    save_full_seq: bool = False,
    save_derot_seq: bool = False,
    adisync: bool = True,
    hwp_inds: tuple[int, int, int, int] = (0, 1, 2, 3),
    extra_pol_angle_corr: float = 0.0,
    plot: bool = True,
) -> PdiResult:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    inst = res.instrument
    dirs = Dirs(cfg.workdir)
    dirs.proc.mkdir(parents=True, exist_ok=True)
    pd = cfg.pdi
    rl = pd.get("rlims", {})
    foip_rlims = tuple(foip_rlims or rl.get("ip", (10, 25)))
    bg_rlims = tuple(bg_rlims or rl.get("bg", (55, 60)))
    unres_pol_rlims = tuple(unres_pol_rlims or rl.get("unres", (0, 6)))
    angoffset = float(angoffset if angoffset is not None else inst.north_offset_deg)
    collapse_type = collapse_type or ("mean" if meanadd else "median")
    nums = cfg.frames("science")
    table = (
        hwp_info
        if isinstance(hwp_info, HwpTable)
        else HwpTable.load(
            Path(hwp_info) if Path(hwp_info).is_absolute() else dirs.root / str(hwp_info),
            inst.pdi.hwp_angles_deg,
        )
    )
    by_num = {n: i for i, n in enumerate(table.filenums)}
    missing = [n for n in nums if n not in by_num]
    if missing:
        raise ValueError(f"frames missing from the HWP table: {missing}")
    if left_suffix is None:
        left_suffix = (
            "leftreg_cal"
            if (dirs.reg / frame_filename(nums[0], "leftreg_cal", prefix)).exists()
            else "leftreg"
        )
    if right_suffix is None:
        right_suffix = left_suffix.replace("left", "right")
    first = read_cube(dirs.reg / frame_filename(nums[0], left_suffix, prefix))
    shape = first.shape
    nlam = first.nlam
    wav = res.wavelengths_um
    center = res.star_center or image_center(shape)
    rad = radial_grid(shape, center)
    phi = polarisation_angle_grid(shape, center)
    terms = MuellerTerms(wav, inst) if (mm_ip_corr or mm_crosstalk_corr) else None
    if first_order_ip_corr and (mm_ip_corr or mm_crosstalk_corr):
        log.warning(
            "first_order_ip_corr is only used without the Mueller-matrix corrections; ignored"
        )
        first_order_ip_corr = False

    # cycles: frame numbers per HWP position
    pos = np.array([table.hwp_pos[by_num[n]] for n in nums])
    cyc = np.array([table.hwp_cycle[by_num[n]] for n in nums])
    ang = np.array([table.hwp_ang[by_num[n]] for n in nums])
    pa_all = {}
    for n in nums:
        h1 = read_cube(dirs.reg / frame_filename(n, left_suffix, prefix)).h1
        pa_all[n] = float(header_get(h1, "PARANG", header_get(h1, "PA", 0.0)))
    cycles = sorted(int(c) for c in np.unique(cyc[cyc >= 0]))
    good_cycles = []
    for c in cycles:
        sel = np.flatnonzero(cyc == c)
        members = {}
        for k, hp in enumerate(hwp_inds):
            m = [nums[i] for i in sel if pos[i] == hp]
            if len(m) != 1:
                break
            members[k] = m[0]
        if len(members) == 4:
            good_cycles.append((c, members))
        else:
            log.warning("cycle %d is incomplete (positions %s); skipped", c, sorted(pos[sel]))
    if not good_cycles:
        raise ValueError("no complete HWP cycles")
    ncyc = len(good_cycles)
    log.info(
        "pdi: %d complete cycles, Mueller corrections ip=%s crosstalk=%s",
        ncyc,
        mm_ip_corr,
        mm_crosstalk_corr,
    )

    q_rows = np.zeros((ncyc, nlam, 3))  # (m_I, m_Q, m_U) of the Q double difference
    u_rows = np.zeros((ncyc, nlam, 3))
    qh = np.full((ncyc, nlam) + shape, np.nan)
    uh = np.full((ncyc, nlam) + shape, np.nan)
    iqh = np.full((ncyc, nlam) + shape, np.nan)
    iuh = np.full((ncyc, nlam) + shape, np.nan)
    q_out_rot = np.zeros(ncyc)
    u_out_rot = np.zeros(ncyc)
    q_hdr = []
    u_hdr = []
    unres_q = np.full((nlam, ncyc), np.nan)
    unres_u = np.full((nlam, ncyc), np.nan)
    cycle_tags = []
    bg_ann = (rad >= bg_rlims[0]) & (rad <= bg_rlims[1])
    foip_ann = (rad >= foip_rlims[0]) & (rad <= foip_rlims[1])
    unres_ann = (rad >= unres_pol_rlims[0]) & (rad <= unres_pol_rlims[1])
    telluric = res.telluric_channels
    first_h0 = first_h1 = None

    for ci, (cnum, members) in enumerate(good_cycles):
        cubes = {}
        for k in range(4):
            n = members[k]
            left = read_cube(dirs.reg / frame_filename(n, left_suffix, prefix))
            right = read_cube(dirs.reg / frame_filename(n, right_suffix, prefix))
            ld, rd = left.data.copy(), right.data.copy()
            if maskspots:
                sats = get_sat_positions(left.h1, nlam)
                for il in range(nlam):
                    for sx, sy in sats[il]:
                        if np.isfinite(sx):
                            m = radial_grid(shape, (sx, sy)) <= 3.0
                            ld[il][m] = np.nan
                            rd[il][m] = np.nan
            cubes[k] = (ld, rd, left.h0, left.h1, n)
        (lqp, rqp, h0_qp, h1_qp, n_qp), (lup, rup, h0_up, h1_up, n_up) = cubes[0], cubes[1]
        (lqm, rqm, h0_qm, h1_qm, n_qm), (lum, rum, h0_um, h1_um, n_um) = cubes[2], cubes[3]
        iq = 0.5 * ((lqp + rqp) + (lqm + rqm))
        iu = 0.5 * ((lup + rup) + (lum + rum))
        pq = 0.5 * ((rqp - lqp) - (rqm - lqm))
        pu = 0.5 * ((rup - lup) - (rum - lum))
        if terms is not None:
            raw_dir = dirs.raw if hasattr(dirs, "raw") else dirs.root / "data" / "raw"
            kw = {"adisync": adisync, "raw_dir": raw_dir}
            xqp = mueller_rows_per_channel(h1_qp, ang[nums.index(n_qp)], terms, h0=h0_qp, **kw)
            xqm = mueller_rows_per_channel(h1_qm, ang[nums.index(n_qm)], terms, h0=h0_qm, **kw)
            xup = mueller_rows_per_channel(h1_up, ang[nums.index(n_up)], terms, h0=h0_up, **kw)
            xum = mueller_rows_per_channel(h1_um, ang[nums.index(n_um)], terms, h0=h0_um, **kw)
            q_rows[ci] = 0.5 * (xqp - xqm)[:, :3]
            u_rows[ci] = 0.5 * (xup - xum)[:, :3]
        # the reference frame of the cycle is the U+ frame (its WCS is used); PA offsets relative to it
        qpa = average_angle([pa_all[n_qp], pa_all[n_qm]])
        upa = average_angle([pa_all[n_up], pa_all[n_um]])
        q_dpa = float(qpa - pa_all[n_up])
        u_dpa = float(upa - pa_all[n_up])
        tag = f"{n_qp:04d},{n_qm:04d}|{n_up:04d},{n_um:04d}"
        cycle_tags.append((cnum, tag))
        # derotate (north-up) each product with its own PA offset
        hq0, hq1 = h0_up.copy(), h1_up.copy()
        qcube, qrot = northup(pq, hq1, angoffset + q_dpa, center=center)
        hu0, hu1 = h0_up.copy(), h1_up.copy()
        ucube, urot = northup(pu, hu1, angoffset + u_dpa, center=center)
        hiq1, hiu1 = h1_up.copy(), h1_up.copy()
        iqcube, _ = northup(iq, hiq1, angoffset + q_dpa, center=center)
        iucube, _ = northup(iu, hiu1, angoffset + u_dpa, center=center)
        q_out_rot[ci], u_out_rot[ci] = qrot, urot
        if save_full_seq:
            for name, data, h0, h1 in (
                ("iq", iq, h0_up, h1_up),
                ("iu", iu, h0_up, h1_up),
                ("q", pq, h0_up, h1_up),
                ("u", pu, h0_up, h1_up),
            ):
                _write(dirs.proc / f"{outfile}_c{cnum:04d}_{name}.fits", data, h0, h1)
        if save_derot_seq:
            _write(dirs.proc / f"{outfile}_c{cnum:04d}_q_derot_uncorr.fits", qcube, hq0, hq1)
            _write(dirs.proc / f"{outfile}_c{cnum:04d}_u_derot_uncorr.fits", ucube, hu0, hu1)
        if first_order_ip_corr:
            for il in range(nlam):
                with np.errstate(divide="ignore", invalid="ignore"):
                    cq = np.nanmedian((qcube[il] / iqcube[il])[foip_ann])
                    cu = np.nanmedian((ucube[il] / iucube[il])[foip_ann])
                qcube[il] -= iqcube[il] * cq
                ucube[il] -= iucube[il] * cu
        if mm_ip_corr and terms is not None:
            qcube -= iqcube * q_rows[ci, :, 0][:, None, None]
            ucube -= iucube * u_rows[ci, :, 0][:, None, None]
        if mm_crosstalk_corr and terms is not None:
            for il in range(nlam):
                a = np.array(
                    [[q_rows[ci, il, 1], q_rows[ci, il, 2]], [u_rows[ci, il, 1], u_rows[ci, il, 2]]]
                )
                try:
                    ainv = np.linalg.inv(a)
                except np.linalg.LinAlgError:
                    log.warning("cycle %d channel %d: singular crosstalk matrix, skipped", cnum, il)
                    continue
                b = np.stack([qcube[il].ravel(), ucube[il].ravel()])
                x = ainv @ b
                qcube[il] = x[0].reshape(shape)
                ucube[il] = x[1].reshape(shape)
        # without HWP/field synchronisation the polarisation axes rotate with the field:
        # undo the north-up rotation (IDL ``outputrot`` = −ROTANG); synchronised data only
        # need the true-north offset
        pol_corr = average_angle([-qrot, -urot]) if not adisync else angoffset
        pol_corr += extra_pol_angle_corr
        if pol_corr != 0:
            qcube, ucube = adjust_pol_axis(qcube, ucube, pol_corr, phi)
        if bg_sub:
            for il in range(nlam):
                qcube[il] -= np.nanmedian(qcube[il][bg_ann])
                ucube[il] -= np.nanmedian(ucube[il][bg_ann])
                iqcube[il] -= np.nanmedian(iqcube[il][bg_ann])
                iucube[il] -= np.nanmedian(iucube[il][bg_ann])
        if unres_pol_sub:
            qrel = np.full(nlam, np.nan)
            urel = np.full(nlam, np.nan)
            for il in range(nlam):
                okq = unres_ann & np.isfinite(qcube[il]) & np.isfinite(iqcube[il])
                oku = unres_ann & np.isfinite(ucube[il]) & np.isfinite(iucube[il])
                if okq.any():
                    qrel[il] = np.nansum(qcube[il][okq]) / np.nansum(iqcube[il][okq])
                if oku.any():
                    urel[il] = np.nansum(ucube[il][oku]) / np.nansum(iucube[il][oku])
            for arr in (qrel, urel):
                fitmask = np.isfinite(arr)
                fitmask[telluric] = False
                bad = ~np.isfinite(arr)
                if bad.any() and fitmask.sum() >= 3:
                    cf = robust_poly_fit(wav[fitmask], arr[fitmask], 2, numit=15)
                    arr[bad] = cf[0] + cf[1] * wav[bad] + cf[2] * wav[bad] ** 2
            unres_q[:, ci], unres_u[:, ci] = qrel, urel
        qh[ci], uh[ci], iqh[ci], iuh[ci] = qcube, ucube, iqcube, iucube
        q_hdr.append((hq0, hq1))
        u_hdr.append((hu0, hu1))
        if ci == 0:
            first_h0, first_h1 = hq0.copy(), hq1.copy()
        if save_derot_seq:
            for name, data, h0, h1 in (
                ("q", qcube, hq0, hq1),
                ("u", ucube, hu0, hu1),
                ("iq", iqcube, hq0, hiq1),
                ("iu", iucube, hu0, hiu1),
            ):
                _write(dirs.proc / f"{outfile}_c{cnum:04d}_{name}_derot.fits", data, h0, h1)

    if unres_pol_sub:
        if avg_unres_pol_coeffs:
            unres_q[:] = np.nanmedian(unres_q, axis=1)[:, None]
            unres_u[:] = np.nanmedian(unres_u, axis=1)[:, None]
        qh -= iqh * unres_q.T[:, :, None, None]
        uh -= iuh * unres_u.T[:, :, None, None]
        np.savetxt(
            dirs.proc / f"{outfile}_unres_pol.txt",
            np.column_stack([wav, unres_q, unres_u]),
            fmt="%.6e",
            header="wavelength_um  q/I per cycle ...  u/I per cycle ...",
        )
        if plot:
            _plot_unres(wav, unres_q, unres_u, dirs.proc / f"{outfile}_unres_pol.png")

    # headers of the combined products
    h0_main = first_h0
    h1_main = first_h1.copy()
    for cnum, tag in cycle_tags:
        h1_main[f"CYC_{cnum:04d}"] = (tag, "frames (Q+,Q-|U+,U-) of this cycle")
    h1_main["NCYCLES"] = (ncyc, "number of HWP cycles combined")
    h1_main["COLLAPSE"] = (collapse_type, "cycle/channel collapse")
    h1_main["ANGOFF"] = (angoffset, "north offset applied [deg]")
    h1_main["MM_IP"] = (int(mm_ip_corr), "Mueller-matrix IP correction")
    h1_main["MM_XTALK"] = (int(mm_crosstalk_corr), "Mueller-matrix crosstalk correction")
    h1_main["FOIPCORR"] = (int(first_order_ip_corr), "first-order IP correction")
    h1_main["BGSUB"] = (int(bg_sub), "background annulus subtraction")
    h1_main["UNRESSUB"] = (int(unres_pol_sub), "unresolved stellar polarisation subtraction")
    h1_main["ADISYNC"] = (int(adisync), "HWP synchronised with field rotation")
    h1_main["HWPINDS"] = (",".join(str(i) for i in hwp_inds), "HWP positions used as Q+,U+,Q-,U-")
    if first_order_ip_corr:
        h1_main["FOIPRLIM"] = (f"{foip_rlims[0]:g},{foip_rlims[1]:g}", "radii [px]")
    if bg_sub:
        h1_main["BG_RLIMS"] = (f"{bg_rlims[0]:g},{bg_rlims[1]:g}", "radii [px]")
    if unres_pol_sub:
        h1_main["UNRRLIMS"] = (f"{unres_pol_rlims[0]:g},{unres_pol_rlims[1]:g}", "radii [px]")
    h1_main["EXPOLANG"] = (extra_pol_angle_corr, "extra polarisation angle correction [deg]")

    files: dict[str, Path] = {}
    if collapse_cycles:
        for ci, (cnum, _) in enumerate(good_cycles):
            qc = _collapse(qh[ci], collapse_type)
            uc = _collapse(uh[ci], collapse_type)
            for name, data in (("q", qc), ("u", uc), ("pi", np.hypot(qc, uc))):
                files[f"c{cnum:04d}_{name}_collapsed"] = _write(
                    dirs.proc / f"{outfile}_c{cnum:04d}_{name}_collapsed.fits",
                    data,
                    h0_main,
                    h1_main,
                )
    qcube = combine_frames(qh, "mean" if meanadd else "median", axis=0)
    ucube = combine_frames(uh, "mean" if meanadd else "median", axis=0)
    pi_cube = np.hypot(qcube, ucube)
    aolp_cube = 0.5 * np.degrees(np.arctan2(ucube, qcube))
    qphi, uphi = azimuthal_stokes(qcube, ucube, phi)
    for name, data in (
        ("q", qcube),
        ("u", ucube),
        ("pi", pi_cube),
        ("aolp", aolp_cube),
        ("qphi", qphi),
        ("uphi", uphi),
    ):
        files[name] = _write(dirs.proc / f"{outfile}_{name}.fits", data, h0_main, h1_main)
    qcol = _collapse(qcube, collapse_type)
    ucol = _collapse(ucube, collapse_type)
    qphicol, uphicol = azimuthal_stokes(qcol, ucol, phi)
    for name, data in (
        ("q", qcol),
        ("u", ucol),
        ("pi", np.hypot(qcol, ucol)),
        ("aolp", 0.5 * np.degrees(np.arctan2(ucol, qcol))),
        ("qphi", qphicol),
        ("uphi", uphicol),
    ):
        files[f"{name}_collapsed"] = _write(
            dirs.proc / f"{outfile}_{name}_collapsed.fits", data, h0_main, h1_main
        )
    coverage = np.mean(np.isfinite(qh + uh).all(axis=1), axis=0)
    files["coverage"] = _write(dirs.proc / "coverage_map.fits", coverage, h0_main, h1_main)
    log.info("pdi: wrote %s_{q,u,pi,aolp,qphi,uphi}[_collapsed].fits to %s", outfile, dirs.proc)
    return PdiResult(
        qcube,
        ucube,
        qphi,
        uphi,
        dirs.proc,
        ncyc,
        files,
        unres_q if unres_pol_sub else None,
        unres_u if unres_pol_sub else None,
    )


def _plot_unres(wav, uq, uu, path: Path) -> None:
    from charis_dpp.plotting import finish, plt

    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    pi = np.hypot(uq, uu)
    aolp = 0.5 * np.degrees(np.arctan2(uu, uq))
    for ax, data, lab in zip(
        axes.ravel(),
        (100 * uq, 100 * uu, 100 * pi, aolp),
        ("% q", "% u", "% pol", "AOLP [deg]"),
        strict=True,
    ):
        for ci in range(data.shape[1]):
            ax.plot(wav, data[:, ci], lw=1, label=f"cycle {ci}")
        ax.set_xlabel("wavelength [um]")
        ax.set_ylabel(lab)
    axes[0, 0].legend(fontsize=6, frameon=False, ncol=2)
    fig.suptitle("unresolved stellar polarisation")
    finish(fig, path)
