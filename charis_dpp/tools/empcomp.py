"""Empirical comparison of a CHARIS spectrum with brown-dwarf spectral libraries
(``tools/empcomp/*.pro``): χ² against the Montreal (Gagné+2015) and Bonnefoy
libraries, Allers & Liu (2013)-style spectral-type/gravity indices, and a
bolometric luminosity estimate.

Fixes vs. IDL: spectral types are parsed with a regex (#249); the indices
follow their wavelength definitions (#250); library entries are matched by
exact file stem (#251); each library reads its own metadata (#252); the
χ² distribution plot has no hand-edited points (#253-#256); the luminosity
helper returns its result (#257).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.astro.units import flam_to
from charis_dpp.config.paths import models_dir
from charis_dpp.io.tables import read_spectrum
from charis_dpp.planetspec.convert import resample_to_charis
from charis_dpp.specphotcal.spec_covar import covariance_matrix
from charis_dpp.tools.atmosfit import good_channels, read_covar_params, read_planet_spectrum
from charis_dpp.tools.chisq import chisq, chisq_covar

log = logging.getLogger(__name__)

__all__ = [
    "spectral_type_number",
    "gravity_class",
    "spt_gravity_indices",
    "LibraryFit",
    "fit_empirical_library",
    "luminosity_from_k",
    "plot_chisq_distribution",
]

_LETTER = {"G": -20.0, "K": -10.0, "M": 0.0, "L": 10.0, "T": 20.0, "Y": 30.0}


def spectral_type_number(spt: str) -> float:
    """``'M8.5'`` → 8.5, ``'L1 low'`` → 11, ``'<T3:'`` → 23, ``'K7'`` → −3; NaN if unparsable."""
    m = re.search(r"([GKMLTY])\s*(\d+(?:\.\d+)?)", str(spt).upper())
    if not m:
        return np.nan
    return _LETTER[m.group(1)] + float(m.group(2))


def gravity_class(spt: str) -> int:
    """0 field, 1 intermediate (β), 2 low (γ), 3 very low (δ) from a Montreal 'spt' string."""
    s = str(spt).lower()
    if "very" in s or "δ" in s or "delta" in s:
        return 3
    if "low" in s or "γ" in s or "gamma" in s:
        return 2
    if "int" in s or "β" in s or "beta" in s:
        return 1
    return 0


def _mean_in(wav, flux, lo, hi) -> float:
    m = (wav >= lo) & (wav <= hi)
    if not m.any():
        return float(np.interp(0.5 * (lo + hi), wav, flux))
    return float(np.nanmean(flux[m]))


def spt_gravity_indices(
    wav_um: np.ndarray, fnu: np.ndarray | None = None, flam: np.ndarray | None = None
) -> dict[str, float]:
    """Allers & Liu (2013)-style indices and the polynomial spectral-type estimates
    (``charis_calcsptgravind.pro``), evaluated by wavelength window."""
    wav = np.asarray(wav_um, dtype=float)
    if flam is None:
        flam = np.asarray(fnu, dtype=float) / wav**2
    f = np.asarray(flam, dtype=float)

    def at(l):
        return float(np.interp(l, wav, f))

    h2o = _mean_in(wav, f, 1.522, 1.575) / _mean_in(wav, f, 1.471, 1.522)
    h2od = _mean_in(wav, f, 1.931, 1.999) / at(2.067)
    h2o1 = at(1.329) / _mean_in(wav, f, 1.284, 1.328)
    h2o2 = at(2.067) / at(2.139)
    h2k = _mean_in(wav, f, 2.14, 2.213) / _mean_in(wav, f, 2.213, 2.29)
    lam1, lamc, lam2 = 1.47, 1.56, 1.67
    f1, fc, f2 = at(lam1), at(lamc), at(lam2)
    hcont = ((lamc - lam1) * f2 / (lam2 - lam1) + (lam2 - lamc) * f1 / (lam2 - lam1)) / fc

    def poly(x, c):
        return c[0] + c[1] * x + c[2] * x**2 + c[3] * x**3

    return {
        "H2O": h2o,
        "H2O-D": h2od,
        "H2O-1": h2o1,
        "H2O-2": h2o2,
        "H2K": h2k,
        "H-cont": hcont,
        "spt_H2O": poly(h2o, (-83.5437, 169.388, -104.424, 24.0476)),
        "spt_H2O-D": poly(h2od, (79.4477, -202.245, 229.884, -97.230)),
        "spt_H2O-1": poly(h2o1, (12.1927, 39.3513, -80.7404, 28.5982)),
        "spt_H2O-2": poly(h2o2, (10.8822, 55.4580, -97.8144, 37.5013)),
    }


@dataclass
class LibraryFit:
    name: str
    spt: str
    spt_number: float
    gravity: int
    chisq: float
    ndof: int
    alpha: float
    indices: dict[str, float]

    @property
    def reduced(self) -> float:
        return self.chisq / max(self.ndof, 1)


def _montreal_meta(libdir: Path) -> dict[str, tuple[str, str]]:
    """file stem → (object name, spectral type) from ``montreal_library_full.tsv``."""
    tsv = libdir / "montreal_library_full.tsv"
    out = {}
    if not tsv.exists():
        return out
    for line in tsv.read_text(errors="replace").splitlines():
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        stems = [c for c in cols if re.match(r"^[A-Za-z0-9+\-_.]+_NIR_", c.strip())]
        stem = stems[-1].strip() if stems else None
        if stem:
            out[stem] = (cols[0].strip(), cols[5].strip())
    return out


def fit_empirical_library(
    spectrum: str | Path,
    *,
    library: str = "montreal",
    library_dir: str | Path | None = None,
    covar: str | Path | None = None,
    band: str | None = None,
    outdir: str | Path | None = None,
    plot: bool = False,
) -> list[LibraryFit]:
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
    if library_dir is None:
        library_dir = (
            models_dir()
            / "planetmodels"
            / ("Montreal" if library == "montreal" else "bonn_speclib")
        )
    libdir = Path(library_dir)
    files = sorted(libdir.glob("*.fits"))
    if not files:
        raise FileNotFoundError(f"no spectra in {libdir}")
    meta = _montreal_meta(libdir) if library == "montreal" else {}
    outdir = Path(outdir or f"{library}_library")
    outdir.mkdir(parents=True, exist_ok=True)
    results = []
    for p in files:
        stem = p.name.split(".")[0]
        try:
            sp = read_spectrum(p)
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: cannot read (%s)", p.name, exc)
            continue
        e_in = sp.error if sp.error is not None else 1e-3 * sp.flux
        pred, perr = resample_to_charis(sp.wavelength_um, sp.flux, wav, e_in)
        pred_mjy = flam_to(wav, pred, "mJy")
        perr_mjy = flam_to(wav, perr, "mJy") if perr is not None else None
        if psi is None:
            r = chisq(flux, err, pred_mjy, err_model=perr_mjy, exclude=exclude)
        else:
            r = chisq_covar(flux, err, pred_mjy, psi, err_model=perr_mjy, exclude=exclude)
        name, spt = meta.get(stem, (stem, ""))
        idx = spt_gravity_indices(wav, fnu=pred_mjy)
        results.append(
            LibraryFit(
                name,
                spt,
                spectral_type_number(spt),
                gravity_class(spt),
                r.chisq,
                r.ndof,
                r.alpha,
                idx,
            )
        )
    results.sort(key=lambda m: m.chisq)
    hdr = "# name spt spt_num grav chisq ndof chisq_red alpha H2O H2O-D H2O-1 H2O-2 H2K H-cont"
    lines = [hdr]
    for m in results:
        i = m.indices
        lines.append(
            f"{m.name.replace(' ', '_')} {m.spt.replace(' ', '_') or '-'} {m.spt_number:.1f} {m.gravity} "
            f"{m.chisq:.3f} {m.ndof} {m.reduced:.3f} {m.alpha:.4e} "
            f"{i['H2O']:.4f} {i['H2O-D']:.4f} {i['H2O-1']:.4f} {i['H2O-2']:.4f} {i['H2K']:.4f} {i['H-cont']:.4f}"
        )
    (outdir / "fitoutcome.txt").write_text("\n".join(lines) + "\n")
    target_idx = spt_gravity_indices(wav, fnu=flux)
    (outdir / "target_indices.txt").write_text(
        "\n".join(f"{k} {v:.5f}" for k, v in target_idx.items()) + "\n"
    )
    if plot and results:
        plot_chisq_distribution(
            outdir / "fitoutcome.txt", outdir / "chisqdist.png", target_indices=target_idx
        )
    log.info(
        "%s library: %d spectra, best %s (chisq_red %.2f)",
        library,
        len(results),
        results[0].name if results else "-",
        results[0].reduced if results else np.nan,
    )
    return results


def luminosity_from_k(
    k_mag: float,
    distance_pc: float,
    spt_number: float | None = None,
    *,
    bc_k: float | None = None,
    m_bol: float | None = None,
) -> tuple[float, float]:
    """``(log L/Lsun, M_bol)`` from a K magnitude with the Liu+2010 BC_K polynomial in the
    numerical spectral type (M0 = 0, L0 = 10), or an explicit ``bc_k``/``m_bol``."""
    if bc_k is None:
        if spt_number is None and m_bol is None:
            raise ValueError("give spt_number, bc_k or m_bol")
        s = spt_number if spt_number is not None else 0.0
        bc_k = 3.9257 - 3.8338e-1 * s + 5.3597e-2 * s**2 - 2.655e-3 * s**3 + 4.0859e-5 * s**4
    if m_bol is None:
        m_bol = bc_k + k_mag - 5 * np.log10(distance_pc / 10.0)
    return float(-0.4 * (m_bol - 4.74)), float(m_bol)


def plot_chisq_distribution(
    fitoutcome: str | Path,
    savefig: str | Path,
    *,
    target_indices: dict | None = None,
    limit: float = 1.8,
) -> None:
    from charis_dpp.plotting import finish, plt

    rows = [
        ln.split()
        for ln in Path(fitoutcome).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    spt = np.array([float(r[2]) for r in rows])
    grav = np.array([int(r[3]) for r in rows])
    chi = np.array([float(r[6]) for r in rows])
    h2k = np.array([float(r[12]) for r in rows])
    hcont = np.array([float(r[13]) for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    cols = {0: "0.6", 1: "C0", 2: "C2", 3: "C1"}
    labels = {0: "field", 1: "β / intermediate", 2: "γ / low", 3: "δ / very low"}
    for g in range(4):
        m = grav == g
        if not m.any():
            continue
        axes[0].scatter(spt[m], chi[m], s=18, color=cols[g], label=labels[g])
        axes[1].scatter(spt[m], hcont[m], s=18, color=cols[g])
        axes[2].scatter(spt[m], h2k[m], s=18, color=cols[g])
    axes[0].set_yscale("log")
    axes[0].axhline(limit, color="k", lw=0.8, ls="--")
    axes[0].set_ylabel("reduced chi^2")
    axes[1].set_ylabel("H-cont")
    axes[2].set_ylabel("H2K")
    if target_indices:
        axes[1].axhline(target_indices["H-cont"], color="C3", lw=1.2, label="target")
        axes[2].axhline(target_indices["H2K"], color="C3", lw=1.2)
    ticks = [-10, -5, 0, 5, 10, 15, 20, 25]
    for ax in axes:
        ax.set_xticks(ticks)
        ax.set_xticklabels(["K0", "K5", "M0", "M5", "L0", "L5", "T0", "T5"])
        ax.set_xlabel("spectral type")
    axes[0].legend(frameon=False, fontsize=8)
    finish(fig, savefig)
