# Python workflow

Commands are `charis <step> <dataset.toml> [options]`; `--help` lists the options (same
names as the IDL keywords). Each step is also a function, e.g.
`from charis_dpp.reg import register_cube; register_cube("HD1160_lowres.toml")`.
Products go under `reduc/{prep,reg,rsub,proc}/`.

## Configuration

* `<Target>_<band>.toml` in the working directory: target, observation (frames, astrogrid,
  modulation), `[reduction]` defaults for every step, `[pdi]`, and an optional `[instrument]`
  table to override a constant for this dataset (e.g. `pixscale_arcsec = 0.0162`).
* `charis_dpp/data/instrument.toml`: instrument constants by date (pixel scale, north offset,
  sat-spot contrasts, NBS parity, Mueller-matrix parameters, M3 recoat...). Picked by the
  observation's MJD; `charis show-config` prints what will be used.
* Command-line options override both for one run.

`charis migrate-info X.info` converts an old parameter file. Without network access,
`newobs` cannot query SIMBAD: fill in `[target]` by hand before `specphot-cal`.

## Total-intensity reduction

| step | command | notes |
|---|---|---|
| 1 | `newobs DATE OBJECT BAND` | creates the `.toml` |
| 2 | `imprep` | headers, science/sky classification (`obslog.txt`) |
| 3 | `subtract-sky [--scalesky]` | |
| 4 | `register [--method sats\|psf\|calsat\|xcorr]` | writes `reduc.log`. `xcorr` needs no satellite spots: frames are aligned to the first one; its centre comes from `--xcorr-center symmetry\|centroid\|given` |
| 5 | `specphot-cal [--calmethod N] [--unsat]`, `make-emp-psf` | flux calibration, PSF model |
| 6 | `imrsub --prad` | spatial filtering |
| 7 | `adialoci`, `sdialoci`, `adiklip`, `rdiklip`, `ssdi` | PSF subtraction → `reduc/proc/final.fits` |
| 8 | `detect`, `fwdmod-planet`, `extract-spectrum`, `attenmap`, `final-contrast` | planets |
| 8' | `disk-fwdmod`, `negspec-disk`, `recalc-disk-chisq` | disks |
| 9 | `atmosfit`, `empcomp`, `spec-covar`, `prop-mot`, `bbphot` | analysis |

`autoreduce [--stop STEP]` runs 1–8.

## PDI

```
charis imprep X.toml --pdi
charis pdi-flatfield X.toml --brightflatdir ../pdi_flats/
charis pdi-subtract-sky X.toml --scalesky
charis pdi-register X.toml            # --method xcorr without satellite spots
charis pdi-specphot-cal X.toml --modamp 25 --ap-factor 3
charis pdi-hwp-match X.toml           # VAMPIRES log, hwp_angles csv, or RET-ANG1 headers
charis pdi X.toml --meanadd [--bg-sub --unres-pol-sub]
```

`pdi-hwp-match` prints one line per frame and warns about frames without an angle, HWP
moves during an exposure, and logs from another night (`--time-ref start|mid|end`
chooses the exposure reference time; default `mid`). Products:
`reduc/proc/<out>_{q,u,pi,aolp,qphi,uphi}[_collapsed].fits`.

## Differences from IDL worth knowing

* arrays are `(channel, y, x)`; header cards keep the IDL `x y` order.
* `ROTANG` is the rotation actually applied for north-up.
* no interactive prompts: a missing input raises an error naming the expected file.
* `BUGS.md` lists every IDL defect and how it is handled.
