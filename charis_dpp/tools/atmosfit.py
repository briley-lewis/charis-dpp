"""Fit a grid of atmosphere models (BT-Settl ``lte*`` files) to a CHARIS spectrum
(``tools/atmosfit/charis_atmosmodelfit.pro``).

Each model (F_λ at the surface, wavelength in Å) is smoothed/resampled to the
CHARIS channels, scaled to ``(R/d)²`` for a 1 R_Jup radius, converted to mJy
and fitted with an analytic scale factor (so ``sqrt(alpha)`` is the radius
in R_Jup); optional spectral covariance and broadband photometry enter the
χ².  Results go to ``<outdir>/fitoutcome.txt`` sorted by χ².

Fixes vs. IDL: the band is derived from the wavelengths (#245); masked
channels come from the telluric ranges (#246); no fixed LUNs, duplicate
reads or ``stop`` (#247).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.astro.units import flam_to
from charis_dpp.astro.zeropoints import bandwidth_um
from charis_dpp.config.paths import btsettl_dir
from charis_dpp.io.wavelengths import band_channels, telluric_channels
from charis_dpp.planetspec.convert import resample_to_charis, synthetic_photometry
from charis_dpp.specphotcal.spec_covar import covariance_matrix
from charis_dpp.tools.chisq import ChiSq, chisq, chisq_covar

log = logging.getLogger(__name__)

__all__ = [
    "ModelFit",
    "fit_model_grid",
    "read_planet_spectrum",
    "read_covar_params",
    "good_channels",
]

R_JUP_PC = 7.1492e7 / 3.0856776e16  # Jupiter radius in parsec


@dataclass
class ModelFit:
    name: str
    teff: float
    logg: float
    metallicity: float
    chisq: float
    ndof: int
    alpha: float
    model_mjy: np.ndarray
    phot_mjy: np.ndarray | None = None

    @property
    def reduced(self) -> float:
        return self.chisq / max(self.ndof, 1)

    @property
    def radius_rjup(self) -> float:
        return float(np.sqrt(self.alpha))


def read_planet_spectrum(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``wavelength_um flux err [snr]`` columns (the ``extract-spectrum`` output)."""
    a = np.loadtxt(path, ndmin=2)
    return a[:, 0], a[:, 1], a[:, 2]


def read_covar_params(path: str | Path) -> dict[str, float]:
    v = np.loadtxt(path, ndmin=2)[0]
    keys = ("a_rho", "a_lam", "a_del", "sig_rho", "sig_lam", "rho")
    return dict(zip(keys, (float(x) for x in v[: len(keys)]), strict=False))


def good_channels(wav: np.ndarray, band: str | None = None) -> np.ndarray:
    """Channels used in the fit: no telluric channels, optionally one sub-band (j/h/k/hk/jhk)."""
    mask = np.ones(wav.size, dtype=bool)
    mask[telluric_channels(wav)] = False
    if band:
        bands = band_channels(wav)
        sel = np.zeros(wav.size, dtype=bool)
        for b in band.upper():
            if b in bands:
                sel[bands[b]] = True
        mask &= sel
    return np.flatnonzero(mask)


def _parse_btsettl_name(name: str) -> tuple[float, float, float]:
    m = re.match(r"lte(\d{3})-?(\d\.\d)([+-]\d\.\d)?", name)
    if not m:
        return np.nan, np.nan, np.nan
    return float(m.group(1)) * 100.0, float(m.group(2)), float(m.group(3) or 0.0)


def _read_model(path: Path) -> tuple[np.ndarray, np.ndarray]:
    a = np.loadtxt(path, ndmin=2, usecols=(0, 1))
    wav = a[:, 0] * 1e-4  # Å → µm
    f = a[:, 1]
    keep = (wav > 0.9) & (wav < 6.0)
    order = np.argsort(wav[keep])
    return wav[keep][order], f[keep][order]


def fit_model_grid(
    spectrum: str | Path,
    *,
    grid_dir: str | Path | None = None,
    models: list[Path] | None = None,
    covar: str | Path | None = None,
    phot: str | Path | None = None,
    band: str | None = None,
    distance_pc: float = 10.0,
    outdir: str | Path = "btsettl_fit",
    plot: bool = False,
) -> list[ModelFit]:
    wav, flux, err = read_planet_spectrum(spectrum)
    good = good_channels(wav, band)
    exclude = np.setdiff1d(np.arange(wav.size), good)
    psi = None
    if covar is not None:
        p = read_covar_params(covar)
        lr = (wav[:, None] - wav[None, :]) / np.median(wav)
        psi = covariance_matrix(
            lr, p["rho"], p["a_rho"], p["a_lam"], p["a_del"], p["sig_rho"], p["sig_lam"]
        )
    phot_names, phot_wav, phot_flux, phot_err = [], None, None, None
    if phot is not None:
        rows = [
            ln.split()
            for ln in Path(phot).read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
        phot_names = [r[0] for r in rows]
        phot_wav = np.array([float(r[1]) for r in rows])
        phot_flux = np.array([float(r[2]) for r in rows])
        phot_err = np.array([float(r[3]) for r in rows])
        for n in phot_names:
            bandwidth_um(n)  # validates the filter name early
    paths = (
        [Path(m) for m in models]
        if models
        else sorted(Path(grid_dir or btsettl_dir()).glob("lte*"))
    )
    if not paths:
        raise FileNotFoundError("no lte* model files found")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scale = (R_JUP_PC / distance_pc) ** 2
    results = []
    for p in paths:
        try:
            mw, mf = _read_model(p)
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: cannot read (%s)", p.name, exc)
            continue
        pred, _ = resample_to_charis(mw, mf, wav)
        pred_mjy = flam_to(wav, pred * scale, "mJy")
        pm = None
        if phot is not None:
            pm = np.array([synthetic_photometry(mw, mf * scale, n)[0] for n in phot_names])
        if psi is None:
            r: ChiSq = chisq(
                flux, err, pred_mjy, exclude=exclude, phot=phot_flux, ephot=phot_err, modelphot=pm
            )
        else:
            r = chisq_covar(
                flux,
                err,
                pred_mjy,
                psi,
                exclude=exclude,
                phot=phot_flux,
                ephot=phot_err,
                modelphot=pm,
            )
        teff, logg, met = _parse_btsettl_name(p.name)
        results.append(
            ModelFit(
                p.name,
                teff,
                logg,
                met,
                r.chisq,
                r.ndof,
                r.alpha,
                pred_mjy * r.alpha,
                None if pm is None else pm * r.alpha,
            )
        )
    results.sort(key=lambda m: m.chisq)
    lines = ["# model teff logg [M/H] chisq ndof chisq_red alpha radius_rjup"]
    for m in results:
        lines.append(
            f"{m.name} {m.teff:.0f} {m.logg:.1f} {m.metallicity:+.1f} {m.chisq:.3f} {m.ndof} {m.reduced:.3f} {m.alpha:.4e} {m.radius_rjup:.3f}"
        )
    (outdir / "fitoutcome.txt").write_text("\n".join(lines) + "\n")
    if plot and results:
        from charis_dpp.plotting import finish, plt

        for m in results[:5]:
            fig, ax = plt.subplots(figsize=(6.5, 4))
            ax.errorbar(wav, flux, err, fmt="o", ms=3, color="C0", label="CHARIS")
            ax.plot(wav, m.model_mjy, "-", color="C3", label=f"{m.name}  χ²ν={m.reduced:.2f}")
            if phot is not None:
                ax.errorbar(phot_wav, phot_flux, phot_err, fmt="s", color="C2", label="photometry")
                ax.plot(phot_wav, m.phot_mjy, "x", color="C3")
            ax.set_xlabel("wavelength [um]")
            ax.set_ylabel("flux density [mJy]")
            ax.legend(frameon=False, fontsize=8)
            finish(fig, outdir / f"{m.name}.png")
    log.info(
        "fit %d models; best %s (chisq_red %.2f)",
        len(results),
        results[0].name if results else "-",
        results[0].reduced if results else np.nan,
    )
    return results
