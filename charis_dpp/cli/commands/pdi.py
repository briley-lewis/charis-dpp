"""PDI subcommands."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register("pdi-flatfield", help="Build the PDI master flat and flat-field the prepared cubes")
def _ff(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--brightflatdir", type=Path)
    parser.add_argument("--darkflatdir", type=Path)
    parser.add_argument("--master", default="master_pdiflat.fits")

    def run(args) -> int:
        from charis_dpp.pdi.flatfield import pdi_flatfield

        out = pdi_flatfield(
            args.config,
            brightflatdir=args.brightflatdir,
            darkflatdir=args.darkflatdir,
            master=args.master,
        )
        print(f"flat-fielded {len(out)} cubes")
        return 0

    return run


@register("pdi-subtract-sky", help="Sky subtraction of the (flat-fielded) PDI cubes")
def _sky(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--scalesky", action="store_true")
    parser.add_argument("--scaleindchan", action="store_true")
    parser.add_argument("--skyrad", type=float)
    parser.add_argument("--skylim", type=float)

    def run(args) -> int:
        from charis_dpp.pdi.subtract_sky import pdi_subtract_sky

        out = pdi_subtract_sky(
            args.config,
            scalesky=args.scalesky,
            scaleindchan=args.scaleindchan,
            skyrad=args.skyrad,
            skylim=args.skylim,
        )
        print(f"sky-subtracted {len(out)} cubes")
        return 0

    return run


@register("pdi-register", help="Register the two Wollaston beams and write the beam sum")
def _reg(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--method", default="sats", help="sats, psf, calsat or xcorr")
    parser.add_argument("--polcent-offset", nargs=2, type=float, metavar=("X", "Y"))
    parser.add_argument("--astrogrid", type=float)
    parser.add_argument("--nomask", action="store_true")
    parser.add_argument("--ladder", action="store_true")
    parser.add_argument("--fitgauss", action="store_true")
    parser.add_argument("--medbox", type=int)

    def run(args) -> int:
        from charis_dpp.pdi.register import pdi_register_cube

        kw = {}
        if args.fitgauss:
            kw["fitgauss"] = True
        if args.medbox is not None:
            kw["medbox"] = args.medbox
        pdi_register_cube(
            args.config,
            method=args.method,
            polcent_offset=tuple(args.polcent_offset) if args.polcent_offset else None,
            astrogrid=args.astrogrid,
            nomask=args.nomask,
            ladder=args.ladder,
            **kw,
        )
        print("wrote reduc/reg/n####{left,right}reg.fits and reduc_left/right.log")
        return 0

    return run


@register("pdi-specphot-cal", help="Spectrophotometric calibration of the two beams")
def _cal(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--calmethod", type=int, default=0)
    parser.add_argument("--unsat", action="store_true")
    parser.add_argument("--modamp", type=float)
    parser.add_argument("--ap-factor", type=float)
    parser.add_argument("--starlib")
    parser.add_argument("--av", type=float)
    parser.add_argument("--empspectrum")
    parser.add_argument(
        "--fov-mask",
        nargs=3,
        type=float,
        metavar=("DX", "DY", "ANGLE"),
        help="full widths (px) and angle (deg) of the valid FOV, as IDL fov_mask_pars",
    )

    def run(args) -> int:
        from charis_dpp.pdi.specphot_cal import pdi_specphot_cal

        out = pdi_specphot_cal(
            args.config,
            fov_mask_pars=tuple(args.fov_mask) if args.fov_mask else None,
            calmethod=args.calmethod,
            unsat=args.unsat,
            modamp=args.modamp,
            ap_factor=args.ap_factor,
            starlib=args.starlib,
            av=args.av,
            empspectrum=args.empspectrum,
        )
        print(f"calibrated {len(out)} beam cubes")
        return 0

    return run


@register("pdi-hwp-match", help="Match frames to HWP angles and cycles (log, csv or headers)")
def _hwp(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--source", choices=["log", "angles", "header"])
    parser.add_argument("--hwp-log", type=Path)
    parser.add_argument("--time-ref", choices=["start", "mid", "end"])
    parser.add_argument("--outfile", default="hwp_info.txt")

    def run(args) -> int:
        from charis_dpp.pdi.hwp_match import hwp_match

        t = hwp_match(
            args.config,
            source=args.source,
            hwp_log=args.hwp_log,
            time_ref=args.time_ref,
            outfile=args.outfile,
        )
        print(f"{len(t.names)} frames, {t.ncycles} complete cycles -> {args.outfile}")
        return 0

    return run


@register("pdi", help="Double-difference polarimetry with Mueller-matrix corrections")
def _pdi(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--hwp-info", default="hwp_info.txt")
    parser.add_argument("--outfile", default="pdi")
    parser.add_argument("--meanadd", action="store_true")
    parser.add_argument("--collapse-type", choices=["median", "mean", "rmean"])
    parser.add_argument("--angoffset", type=float)
    parser.add_argument("--maskspots", action="store_true")
    parser.add_argument("--no-mm-ip-corr", action="store_true")
    parser.add_argument("--no-mm-crosstalk-corr", action="store_true")
    parser.add_argument("--first-order-ip-corr", action="store_true")
    parser.add_argument("--foip-rlims", nargs=2, type=float)
    parser.add_argument("--bg-sub", action="store_true")
    parser.add_argument("--bg-rlims", nargs=2, type=float)
    parser.add_argument("--unres-pol-sub", action="store_true")
    parser.add_argument("--unres-pol-rlims", nargs=2, type=float)
    parser.add_argument("--avg-unres-pol-coeffs", action="store_true")
    parser.add_argument("--collapse-cycles", action="store_true")
    parser.add_argument("--save-full-seq", action="store_true")
    parser.add_argument("--save-derot-seq", action="store_true")
    parser.add_argument("--no-adisync", action="store_true")
    parser.add_argument("--hwp-inds", nargs=4, type=int, default=(0, 1, 2, 3))
    parser.add_argument("--extra-pol-angle-corr", type=float, default=0.0)

    def run(args) -> int:
        from charis_dpp.pdi.pdi import pdi

        r = pdi(
            args.config,
            hwp_info=args.hwp_info,
            outfile=args.outfile,
            meanadd=args.meanadd,
            collapse_type=args.collapse_type,
            angoffset=args.angoffset,
            maskspots=args.maskspots,
            mm_ip_corr=not args.no_mm_ip_corr,
            mm_crosstalk_corr=not args.no_mm_crosstalk_corr,
            first_order_ip_corr=args.first_order_ip_corr,
            foip_rlims=tuple(args.foip_rlims) if args.foip_rlims else None,
            bg_sub=args.bg_sub,
            bg_rlims=tuple(args.bg_rlims) if args.bg_rlims else None,
            unres_pol_sub=args.unres_pol_sub,
            unres_pol_rlims=tuple(args.unres_pol_rlims) if args.unres_pol_rlims else None,
            avg_unres_pol_coeffs=args.avg_unres_pol_coeffs,
            collapse_cycles=args.collapse_cycles,
            save_full_seq=args.save_full_seq,
            save_derot_seq=args.save_derot_seq,
            adisync=not args.no_adisync,
            hwp_inds=tuple(args.hwp_inds),
            extra_pol_angle_corr=args.extra_pol_angle_corr,
        )
        print(f"{r.ncycles} cycles -> {r.outdir}/{args.outfile}_*.fits")
        return 0

    return run


@register("pdi-rlims-grid", help="Grid search over the unresolved-polarisation annulus radii")
def _grid(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--rmin", nargs="+", type=float, default=(0, 2, 4, 6, 8, 10))
    parser.add_argument("--rmax", nargs="+", type=float, default=(5, 10, 15, 20, 25, 30))

    def run(args) -> int:
        from charis_dpp.pdi.rlims_gridsearch import rlims_gridsearch

        runs = rlims_gridsearch(args.config, args.rmin, args.rmax)
        print(f"{len(runs)} runs")
        return 0

    return run
