"""``charis specphot-cal`` / ``satspot-cal`` / ``make-emp-psf`` / ``make-synth-psf`` / ``bbphot``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register("specphot-cal", help="Spectrophotometric calibration of the registered cubes")
def _spc(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--calmethod",
        type=int,
        default=0,
        help="0 each cube | 1 median cube | 2 --calcube | 3 ladder attenuation | 4 ladder cubes | 5 copy from another star",
    )
    parser.add_argument(
        "--unsat", action="store_true", help="unsaturated star instead of satellite spots"
    )
    parser.add_argument("--calcube", nargs="+", type=Path)
    parser.add_argument(
        "--datacube", nargs="+", type=Path, help="cubes to calibrate (default: science frames)"
    )
    parser.add_argument("--modamp", type=float, help="spot modulation amplitude [nm]")
    parser.add_argument(
        "--atten-fac", type=float, help="spot/star flux ratio to use instead of the epoch table"
    )
    parser.add_argument("--subskyannulus", action="store_true")
    parser.add_argument(
        "--spot", choices=["left", "right", "upperleft", "lowerleft", "upperright", "lowerright"]
    )
    parser.add_argument("--fitbckgd", action="store_true")
    parser.add_argument("--expandederror", action="store_true")
    parser.add_argument("--filtcal-slice", action="store_true")
    parser.add_argument("--nopradcal", action="store_true")
    parser.add_argument("--meancomb", action="store_true")
    parser.add_argument("--fixradius", type=float)
    parser.add_argument("--ap-factor", type=float)
    parser.add_argument("--filtername", help="normalisation band override (J/H/K)")
    parser.add_argument("--starlib", choices=["kurucz", "pickles", "empirical"])
    parser.add_argument("--empspectrum", type=Path)
    parser.add_argument("--av", type=float)
    parser.add_argument(
        "--fluxunits",
        default="mJy",
        help="mJy | Jy | W/m^2/um | ergs/s/cm^2/A | ergs/s/cm^2/Hz | contrast",
    )
    parser.add_argument("--suffix", default="reg")

    def run(args) -> int:
        from charis_dpp.specphotcal import specphot_cal

        out = specphot_cal(
            args.config,
            calmethod=args.calmethod,
            unsat=args.unsat,
            calcube=args.calcube,
            datacube=args.datacube,
            modamp=args.modamp,
            atten_fac=args.atten_fac,
            subskyannulus=args.subskyannulus,
            spot_select=args.spot,
            fitbckgd=args.fitbckgd,
            expandederror=args.expandederror,
            filtcal_slice=args.filtcal_slice,
            nopradcal=args.nopradcal,
            meancomb=args.meancomb,
            fixradius=args.fixradius,
            ap_factor=args.ap_factor,
            filtername=args.filtername,
            starlib=args.starlib,
            empspectrum=args.empspectrum,
            av=args.av,
            fluxunits=args.fluxunits,
            suffix=args.suffix,
        )
        print(f"calibrated {len(out)} cubes")
        return 0

    return run


@register("satspot-cal", help="Fit the satellite-spot contrast from ladder cubes")
def _ssc(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--calmethod", type=int, default=0)
    parser.add_argument("--siglim", type=float)
    parser.add_argument("--subskyannulus", action="store_true")
    parser.add_argument("--nopradcal", action="store_true")
    parser.add_argument("--meancomb", action="store_true")
    parser.add_argument("--ap-factor", type=float)

    def run(args) -> int:
        from charis_dpp.specphotcal import satspot_cal

        p0, _ = satspot_cal(
            args.config,
            calmethod=args.calmethod,
            siglim=args.siglim,
            subskyannulus=args.subskyannulus,
            nopradcal=args.nopradcal,
            meancomb=args.meancomb,
            ap_factor=args.ap_factor,
        )
        print(f"spot/star contrast at reference wavelength: {p0:.4e}")
        return 0

    return run


@register("make-emp-psf", help="Build the empirical PSF cube from the satellite spots")
def _emp(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--psfsize", type=int)
    parser.add_argument("--fsize", type=int)
    parser.add_argument("--unsat", action="store_true")
    parser.add_argument("--filt", action="store_true")
    parser.add_argument("--nozero", action="store_true")
    parser.add_argument("--nonorm", action="store_true")
    parser.add_argument("--ladder", action="store_true")
    parser.add_argument("--suffix")
    parser.add_argument("--outname")

    def run(args) -> int:
        from charis_dpp.specphotcal import make_empirical_psf

        p = make_empirical_psf(
            args.config,
            psfsize=args.psfsize,
            fsize=args.fsize,
            unsat=args.unsat,
            filt=args.filt,
            nozero=args.nozero,
            nonorm=args.nonorm,
            ladder=args.ladder,
            suffix=args.suffix,
            outname=args.outname,
        )
        print(f"wrote {p}")
        return 0

    return run


@register("make-synth-psf", help="Build a synthetic (Gaussian or flat) PSF cube")
def _syn(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--kind", choices=["gaussian", "flat"], default="gaussian")
    parser.add_argument("--fwhm", type=float, help="pixels (or arcsec if < 1)")
    parser.add_argument("--rmax", type=float)
    parser.add_argument("--psfsize", type=int)
    parser.add_argument("--outname")

    def run(args) -> int:
        from charis_dpp.specphotcal import make_synthetic_psf

        p = make_synthetic_psf(
            args.config,
            kind=args.kind,
            fwhm=args.fwhm,
            rmax=args.rmax,
            psfsize=args.psfsize,
            outname=args.outname,
        )
        print(f"wrote {p}")
        return 0

    return run


@register(
    "bbphot", help="Broadband J/H/Ks photometry of a spectrum file (wavelength[um] flux [err])"
)
def _bb(parser):
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("--units", default="mJy")

    def run(args) -> int:
        from charis_dpp.specphotcal.bbphot import broadband_photometry_file

        for b, r in broadband_photometry_file(args.spectrum, args.units).items():
            print(
                f"{b:3s} {r.flux_mjy:12.4f} +/- {r.err_mjy:10.4f} mJy  at {r.lam_eff_um:.3f} um   mag {r.mag:.3f} +/- {r.mag_err:.3f}"
            )
        return 0

    return run
