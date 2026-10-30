"""Matplotlib helpers.  Every plotting function takes ``savefig`` and never opens a window."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

__all__ = ["plot_spectrum", "plot_image", "plot_contrast", "finish"]


def finish(fig, savefig: str | Path | None, close: bool = True) -> Path | None:
    out = None
    if savefig:
        out = Path(savefig)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight", dpi=150)
    if close:
        plt.close(fig)
    return out


def plot_spectrum(
    wav_um: np.ndarray,
    flux: np.ndarray,
    err: np.ndarray | None = None,
    *,
    label: str | None = None,
    units: str = "mJy",
    overplot: list[tuple[np.ndarray, np.ndarray, str]] | None = None,
    title: str | None = None,
    savefig: str | Path | None = None,
    from_zero: bool = False,
):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    if err is not None:
        ax.errorbar(wav_um, flux, yerr=err, fmt="o-", ms=4, capsize=2, label=label)
    else:
        ax.plot(wav_um, flux, "o-", ms=4, label=label)
    for w, f, lab in overplot or []:
        ax.plot(w, f, "-", lw=1, label=lab)
    ax.set_xlabel(r"Wavelength ($\mu$m)")
    ax.set_ylabel(f"$F_\\nu$ ({units})")
    if from_zero:
        ax.set_ylim(bottom=0)
    if title:
        ax.set_title(title)
    if label or overplot:
        ax.legend(frameon=False)
    return finish(fig, savefig)


def plot_image(
    im: np.ndarray,
    *,
    title: str | None = None,
    vmin=None,
    vmax=None,
    points=None,
    savefig=None,
    cmap="magma",
):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(im, origin="lower", vmin=vmin, vmax=vmax, cmap=cmap)
    if points is not None:
        pts = np.atleast_2d(points)
        ax.plot(pts[:, 0], pts[:, 1], "c+", ms=10)
    if title:
        ax.set_title(title)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    return finish(fig, savefig)


def plot_contrast(
    sep_arcsec: np.ndarray, curves: dict[str, np.ndarray], *, title=None, savefig=None
):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for lab, c in curves.items():
        ax.plot(sep_arcsec, c, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("Separation (arcsec)")
    ax.set_ylabel(r"5$\sigma$ contrast")
    ax.legend(frameon=False, fontsize=8)
    if title:
        ax.set_title(title)
    return finish(fig, savefig)
