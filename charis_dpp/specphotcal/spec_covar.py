"""Spectral covariance of the residual noise (``charis_calc_spec_covar.pro``; Greco & Brandt 2016).

The correlation ψ_ij between channels i and j in an annulus around the
planet separation ``rho`` (in λ/D) is measured and fitted with

    ψ(x) = A_ρ exp(-½ (ρ x / σ_ρ)²) + A_λ exp(-½ (x / σ_λ)²) + A_δ δ_ij,
    x = (λ_i − λ_j) / λ_c,   A_δ = 1 − A_ρ − A_λ.

Fixes vs. IDL: the model function has no undefined variable (#110); the
zero-lag value is ``A_δ + A_ρ + A_λ`` (#111); bins do not overlap and the
spline/average arrays have the same length (#112/#113); sub-band index
sets come from the wavelength table (#114).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.filters import idl_filter_image, radial_profile_2d
from charis_dpp.io.fits import read_cube, write_image
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["CovarResult", "spectral_covariance", "gb_model", "covariance_matrix"]


@dataclass
class CovarResult:
    rho: float
    a_rho: float
    a_lam: float
    a_del: float
    sig_rho: float
    sig_lam: float
    psi: np.ndarray  # (nλ, nλ) measured correlation
    lam_ratio: np.ndarray  # (nλ, nλ) (λi-λj)/λc
    params_file: Path | None = None

    def model(self) -> np.ndarray:
        return covariance_matrix(
            self.lam_ratio, self.rho, self.a_rho, self.a_lam, self.a_del, self.sig_rho, self.sig_lam
        )


def gb_model(x: np.ndarray, p: np.ndarray, rho: float) -> np.ndarray:
    """Greco & Brandt correlation vs ``x = Δλ/λ_c`` for ``p = (A_ρ, A_λ, σ_ρ, σ_λ)``."""
    return p[0] * np.exp(-0.5 * (rho * x / p[2]) ** 2) + p[1] * np.exp(-0.5 * (x / p[3]) ** 2)


def covariance_matrix(lam_ratio, rho, a_rho, a_lam, a_del, sig_rho, sig_lam) -> np.ndarray:
    m = gb_model(lam_ratio, np.array([a_rho, a_lam, sig_rho, sig_lam]), rho)
    return m + a_del * np.eye(lam_ratio.shape[0])


def spectral_covariance(
    cfg: DatasetConfig | str | Path,
    datacube: str | Path = "final.fits",
    *,
    rho: float | None = None,
    width: float | None = None,
    step: float | None = None,
    prad: bool = False,
    filt: bool = False,
    normnoise: bool = False,
    maskcoord: tuple[float, float] | None = None,
    excval: float | None = None,
    maskresspot: bool = False,
    fixpos: tuple[float, float] | None = None,
    plot: bool = True,
) -> CovarResult:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    cv = cfg.reduction.get("covar", {})
    rho = float(rho if rho is not None else cv.get("rho_px", 10.0))
    width = float(width if width is not None else cv.get("width", 0.5))
    step = float(step if step is not None else cv.get("step", 0.125))
    excval = float(excval if excval is not None else cv.get("excval", 3.0))
    path = (
        Path(datacube)
        if Path(datacube).is_absolute() or Path(datacube).exists()
        else dirs.proc / datacube
    )
    c = read_cube(path)
    cube = c.data.copy()
    nlam = c.nlam
    wav = res.wavelengths_um
    fwhm = res.fwhm_px
    med_fwhm = float(np.median(fwhm))
    lam_med = float(np.median(wav))
    center = res.star_center or image_center(c.shape)
    r = radial_grid(c.shape, center)
    rho_px = rho * med_fwhm  # rho in λ/D → pixels
    roi = (r >= (rho - width) * med_fwhm) & (r <= (rho + width) * med_fwhm)
    if maskcoord is not None:
        roi &= radial_grid(c.shape, maskcoord) > excval * med_fwhm
    for il in range(nlam):
        im = cube[il]
        if prad:
            im = im - radial_profile_2d(im, 1.0, 1.0, min(c.shape) // 2, center)
        if filt:
            im = im - idl_filter_image(im, median=int(round(10 * med_fwhm)))
        if normnoise:
            im = im / radial_profile_2d(np.abs(im) / 0.6745, 1.0, 1.0, min(c.shape) // 2, center)
        if maskresspot:
            im = np.where(im > 10.0, np.nan, im)
        cube[il] = im
    psi = np.zeros((nlam, nlam))
    lam_ratio = (wav[:, None] - wav[None, :]) / lam_med
    sl = cube[:, roi]
    for i in range(nlam):
        for j in range(nlam):
            num = np.nanmean(sl[i] * sl[j])
            den = np.sqrt(np.nanmean(sl[i] ** 2) * np.nanmean(sl[j] ** 2))
            psi[i, j] = abs(num / den) if den > 0 else np.nan
    # bin the off-diagonal correlation in rho*Δλ/λc
    x = rho * lam_ratio
    edges = np.arange(np.floor(x.min()), np.ceil(x.max()) + step, step)
    centres = 0.5 * (edges[:-1] + edges[1:])
    avg = np.full(centres.size, np.nan)
    err = np.full(centres.size, np.nan)
    off = ~np.eye(nlam, dtype=bool)
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        sel = off & (x >= lo) & (x < hi) & np.isfinite(psi)
        if sel.sum() >= 2:
            avg[k] = np.mean(psi[sel])
            err[k] = np.std(psi[sel])
    good = np.isfinite(avg) & (np.abs(centres) >= step)
    err_med = np.nanmedian(err[good]) if good.any() else 1.0
    p0 = np.array([0.1, 0.3, 0.38, 0.75])
    lb = np.array([5e-4, 5e-4, 0.01, 0.025])
    ub = np.array([0.9, 0.9, 0.9, 2.9])

    def resid(p):
        return (gb_model(centres[good] / rho, p, rho) - avg[good]) / max(err_med, 1e-3)

    fit = least_squares(resid, p0, bounds=(lb, ub))
    a_rho, a_lam, sig_rho, sig_lam = fit.x
    a_del = 1.0 - a_rho - a_lam
    out = CovarResult(
        rho,
        float(a_rho),
        float(a_lam),
        float(a_del),
        float(sig_rho),
        float(sig_lam),
        psi,
        lam_ratio,
    )
    tag = f"{rho:g}".replace(".", "p")
    write_image(dirs.proc / f"spec_covar{tag}lam_D.fits", psi)
    pf = dirs.proc / "covar_params"
    pf.write_text(f"{a_rho:.6f} {a_lam:.6f} {a_del:.6f} {sig_rho:.6f} {sig_lam:.6f} {rho:.6f}\n")
    out.params_file = pf
    if plot:
        from charis_dpp.plotting import finish, plt

        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.plot(x[off].ravel(), psi[off].ravel(), ".", color="0.7", ms=3, label="pairs")
        ax.errorbar(centres[good], avg[good], err[good], fmt="o", color="C1", label="binned")
        xx = np.linspace(centres[good].min(), centres[good].max(), 300)
        ax.plot(
            xx,
            gb_model(xx / rho, fit.x, rho),
            "-",
            color="C3",
            label=f"fit: Aρ={a_rho:.2f} Aλ={a_lam:.2f} Aδ={a_del:.2f}",
        )
        ax.set_xlabel(r"$\rho\,(\lambda_i-\lambda_j)/\lambda_c$")
        ax.set_ylabel(r"$\psi_{ij}$")
        ax.set_title(f"spectral correlation at {rho:g} λ/D ({rho_px:.1f} px)")
        ax.legend(frameon=False, fontsize=8)
        finish(fig, dirs.proc / f"spec_covar{tag}lam_D.png")
    log.info(
        "spectral covariance at %.1f lambda/D: A_rho %.3f A_lam %.3f A_del %.3f sig_rho %.3f sig_lam %.3f",
        rho,
        a_rho,
        a_lam,
        a_del,
        sig_rho,
        sig_lam,
    )
    return out
