"""``charis fwdmod-planet`` / ``fake-planet`` / ``attenmap`` / ``negspec-planet`` / ``extract-spectrum``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


def _planet_args(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--reducname", default="final.fits", help="reduced cube in reduc/proc")
    parser.add_argument(
        "--targfile", type=Path, help="text file: x y log10(contrast) per planet (north-up frame)"
    )
    parser.add_argument(
        "--planet",
        nargs=3,
        type=float,
        metavar=("X", "Y", "LOGC"),
        action="append",
        help="planet position and log contrast (repeatable)",
    )
    parser.add_argument(
        "--psf", default="empirical", help="empirical | gaussian | path to a PSF cube"
    )
    parser.add_argument(
        "--planet-model", help="model spectrum (file name in models/planetmodels or path)"
    )
    parser.add_argument(
        "--reductype", type=int, default=0, help="0 ADI | 1 SDI | 2 ADI then SDI | 3 SDI then ADI"
    )
    parser.add_argument("--sdi-reducname")
    parser.add_argument("--refpath", type=Path)


def _planets(args):
    from charis_dpp.partialsub.inject import Planet

    if args.planet:
        return [Planet(x, y, c) for x, y, c in args.planet]
    return None


@register(
    "fwdmod-planet", help="Forward-model planets through the PSF subtraction to measure throughput"
)
def _fwd(parser):
    _planet_args(parser)
    parser.add_argument("--outfile", default="res0.fits")

    def run(args) -> int:
        from charis_dpp.partialsub.fwdmod_planet import fwdmod_planet

        thr = fwdmod_planet(
            args.config,
            reducname=args.reducname,
            planets=_planets(args),
            targfile=args.targfile,
            psf=args.psf,
            planet_model=args.planet_model,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
            outfile=args.outfile,
        )
        for j, t in enumerate(thr):
            print(
                f"planet {j}: median throughput {float(__import__('numpy').nanmedian(t)):.3f}  -> reduc/proc/synth_throughput{j}.txt"
            )
        return 0

    return run


@register("fake-planet", help="Inject fake planets into the data and re-reduce")
def _fake(parser):
    _planet_args(parser)
    parser.add_argument("--outfile", default="fake0.fits")

    def run(args) -> int:
        from charis_dpp.partialsub.fakeplanet import fake_planet

        fake_planet(
            args.config,
            reducname=args.reducname,
            planets=_planets(args),
            targfile=args.targfile,
            psf=args.psf,
            planet_model=args.planet_model,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
            outfile=args.outfile,
        )
        print("wrote reduc/proc/fakeplanet_res.txt")
        return 0

    return run


@register("attenmap", help="Attenuation (throughput) map from a grid of forward-modelled planets")
def _att(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--reducname", default="final.fits")
    parser.add_argument("--contrast", type=float, help="log10 contrast of the grid planets")
    parser.add_argument("--ntc", type=int)
    parser.add_argument("--lrmin", type=float)
    parser.add_argument("--lrmax", type=float)
    parser.add_argument("--dlrstep", type=float)
    parser.add_argument("--psf", default="empirical")
    parser.add_argument("--planet-model")
    parser.add_argument("--reductype", type=int, default=0)
    parser.add_argument("--sdi-reducname")
    parser.add_argument("--refpath", type=Path)

    def run(args) -> int:
        from charis_dpp.partialsub.attenmap import attenuation_map

        pl, pc = attenuation_map(
            args.config,
            reducname=args.reducname,
            contrast=args.contrast,
            ntc=args.ntc,
            lrmin=args.lrmin,
            lrmax=args.lrmax,
            dlrstep=args.dlrstep,
            psf=args.psf,
            planet_model=args.planet_model,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
        )
        print(f"wrote {pl} and {pc}")
        return 0

    return run


@register(
    "negspec-planet",
    help="Inject the negative of an extracted spectrum at the planet position and re-reduce",
)
def _neg(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--reducname", default="final.fits")
    parser.add_argument("--coords", nargs=2, type=float, required=True, metavar=("X", "Y"))
    parser.add_argument("--spectrum", default="spectrum.dat")
    parser.add_argument("--psf", default="empirical")
    parser.add_argument("--reductype", type=int, default=0)
    parser.add_argument("--sdi-reducname")
    parser.add_argument("--refpath", type=Path)

    def run(args) -> int:
        from charis_dpp.partialsub.negspec_planet import negative_planet

        negative_planet(
            args.config,
            reducname=args.reducname,
            planet=tuple(args.coords),
            spectrum_file=args.spectrum,
            psf=args.psf,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
        )
        print("wrote reduc/proc/res0_np.fits and negplanet_residual.txt")
        return 0

    return run


@register("extract-spectrum", help="Extract a 1-D spectrum at a position in a reduced cube")
def _ext(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--datacube", default="final.fits")
    parser.add_argument("--coords", nargs=2, type=float, required=True, metavar=("X", "Y"))
    parser.add_argument("--throughput", help="synth_throughput<i>.txt to divide by")
    parser.add_argument("--no-fitcollapse", action="store_true")
    parser.add_argument("--fitslice", action="store_true")
    parser.add_argument("--fitbckgd", action="store_true")
    parser.add_argument("--extended", action="store_true")
    parser.add_argument("--fitrange", nargs=2, type=int, metavar=("FIRST", "LAST"))
    parser.add_argument("--nozerosky", action="store_true")
    parser.add_argument("--noprad", action="store_true")
    parser.add_argument("--filt-slice", action="store_true")
    parser.add_argument("--filtsnr", action="store_true")
    parser.add_argument("--nocalerror", action="store_true")
    parser.add_argument("--bklim", type=float, default=4.0)
    parser.add_argument("--bklir", type=float, default=2.0)
    parser.add_argument("--disk", action="store_true", help="surface brightness per arcsec^2")
    parser.add_argument("--outputspec")

    def run(args) -> int:
        from charis_dpp.specphotcal.extract_spectrum import extract_spectrum

        fr = list(range(args.fitrange[0], args.fitrange[1] + 1)) if args.fitrange else None
        r = extract_spectrum(
            args.config,
            args.datacube,
            tuple(args.coords),
            fitcollapse=not args.no_fitcollapse,
            fitslice=args.fitslice,
            fitbckgd=args.fitbckgd,
            extended=args.extended,
            throughput=args.throughput,
            fitrange=fr,
            nozerosky=args.nozerosky,
            noprad=args.noprad,
            filt_slice=args.filt_slice,
            filtsnr=args.filtsnr,
            nocalerror=args.nocalerror,
            bklim=args.bklim,
            bklir=args.bklir,
            mode="disk" if args.disk else "point",
            outputspec=args.outputspec,
        )
        print(
            f"position ({r.position[0]:.2f}, {r.position[1]:.2f}); collapsed SNR {r.snr_collapsed:.1f}; wrote {r.path}"
        )
        return 0

    return run


def _roi_args(parser):
    parser.add_argument(
        "--roirange",
        nargs=3,
        type=float,
        metavar=("DX", "DY", "ANGLE"),
        help="rectangular evaluation region (size and position angle)",
    )
    parser.add_argument("--roioffset", nargs=2, type=float, default=(0.0, 0.0), metavar=("X", "Y"))
    parser.add_argument("--lrmin", type=float)
    parser.add_argument("--lrmax", type=float)


@register("disk-fwdmod", help="Forward-model disk images through the PSF subtraction")
def _dfm(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("models", nargs="+", type=Path, help="disk model FITS files (north-up)")
    parser.add_argument("--reducname", default="final.fits")
    parser.add_argument("--psf", default="empirical", help="'empirical', 'gaussian' or a PSF cube")
    parser.add_argument("--noconv", action="store_true", help="models are already PSF-convolved")
    parser.add_argument("--absflux", action="store_true", help="do not pre-scale the model")
    parser.add_argument(
        "--synthscale", type=float, help="initial model scale instead of the LS guess"
    )
    parser.add_argument("--simplescale", action="store_true", help="unweighted final scaling")
    parser.add_argument("--manscale", action="store_true", help="no final scaling (scale = 1)")
    parser.add_argument("--snrlim", type=float, help="restrict the ROI to SNR > this")
    parser.add_argument("--maskrad", type=float, default=0.0)
    parser.add_argument(
        "--calcattenfact", action="store_true", help="also per-channel scaling and attenuation"
    )
    parser.add_argument("--reductype", type=int, default=0)
    parser.add_argument("--sdi-reducname")
    parser.add_argument("--refpath", type=Path)
    parser.add_argument("--nsc", type=int)
    _roi_args(parser)

    def run(args) -> int:
        from charis_dpp.partialsub.disk_fwdmod import disk_fwdmod

        results = disk_fwdmod(
            args.config,
            args.models,
            reducname=args.reducname,
            psf=args.psf,
            noconv=args.noconv,
            absflux=args.absflux,
            synthscale=args.synthscale,
            simplescale=args.simplescale,
            manscale=args.manscale,
            roirange=tuple(args.roirange) if args.roirange else None,
            roioffset=tuple(args.roioffset),
            snrlim=args.snrlim,
            lrmin=args.lrmin,
            lrmax=args.lrmax,
            maskrad=args.maskrad,
            calcattenfact=args.calcattenfact,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
            nsc=args.nsc,
        )
        for r in results:
            print(
                f"{r.model:30s} scale {r.scale:.4g}  chi2_red {r.chisq_binned:.3f}  ({r.nbins} bins)"
            )
        return 0

    return run


@register("negspec-disk", help="Subtract a disk model from every frame and re-reduce")
def _nd(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("diskmodel", type=Path, help="model image or cube to subtract (north-up)")
    parser.add_argument("--reducname", default="final.fits")
    parser.add_argument("--reductype", type=int, default=0)
    parser.add_argument("--sdi-reducname")
    parser.add_argument("--refpath", type=Path)
    parser.add_argument("--outfile", default="res0_np.fits")
    parser.add_argument("--nsc", type=int)

    def run(args) -> int:
        from charis_dpp.partialsub.negspec_disk import negative_disk

        p = negative_disk(
            args.config,
            args.diskmodel,
            reducname=args.reducname,
            reductype=args.reductype,
            sdi_reducname=args.sdi_reducname,
            refpath=args.refpath,
            outfile=args.outfile,
            nsc=args.nsc,
        )
        print(f"wrote {p}")
        return 0

    return run


@register("recalc-disk-chisq", help="Recompute the chi^2 of forward-modelled disks with a new ROI")
def _rdc(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--reducname", default="final.fits")
    parser.add_argument("--models", nargs="*", type=Path, help="default: reduc/proc/model_psfsub/*")
    parser.add_argument("--snoise", type=float, default=1.0, help="divide the noise map by this")
    _roi_args(parser)

    def run(args) -> int:
        from charis_dpp.partialsub.recalc_disk_chisq import recalc_disk_chisq

        for name, chi, n in recalc_disk_chisq(
            args.config,
            reducname=args.reducname,
            models=args.models or None,
            roirange=tuple(args.roirange) if args.roirange else None,
            roioffset=tuple(args.roioffset),
            lrmin=args.lrmin if args.lrmin is not None else 5.0,
            lrmax=args.lrmax,
            snoise=args.snoise,
        ):
            print(f"{name:30s} chi2_red {chi:.3f} ({n} bins)")
        return 0

    return run
