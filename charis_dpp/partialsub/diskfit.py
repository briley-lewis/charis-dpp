"""Shared pieces of the disk forward-modelling drivers: region of interest,
model scaling, χ² maps and the model-parameter bookkeeping
(``charis_generate_roi.pro``, ``charis_calc_chisqbin.pro`` and the
~300 lines each disk driver repeated).

Fixes vs. IDL: the ROI rectangle is sliced symmetrically for even sizes and
an empty mask is handled (#281); ``snrlim`` is actually used (#179); the χ²
binning is a block mean, not nearest-neighbour sampling (#279); zero noise
pixels are excluded from the χ² instead of giving Inf (#185).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io.fits import Header
from scipy.signal import fftconvolve

from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.geometry.masks import roi_rect_rotated
from charis_dpp.image.interp import rotate_image
from charis_dpp.partialsub.inject import smear_offsets

log = logging.getLogger(__name__)

__all__ = [
    "MODEL_KEYS",
    "DiskModel",
    "DiskFitResult",
    "load_disk_model",
    "roi_mask",
    "scale_model_ls",
    "chisq_map",
    "convolve_model",
    "model_in_frame",
    "fit_outcomes_line",
]

# GRaTeR-style disk parameters carried in the model FITS header (any may be absent)
MODEL_KEYS = (
    "G",
    "KSI0",
    "ALPHAIN",
    "ALPHAOUT",
    "BETA",
    "R0",
    "XDO",
    "YDO",
    "E",
    "THETA0",
    "ITILT",
)


@dataclass
class DiskModel:
    path: Path
    data: np.ndarray  # (ny, nx) or (nλ, ny, nx), north-up, same geometry as the reduced cube
    params: dict[str, float] = field(default_factory=dict)

    @property
    def name(self) -> str:
        return self.path.name.split(".")[0]

    @property
    def is_cube(self) -> bool:
        return self.data.ndim == 3

    def channel(self, il: int) -> np.ndarray:
        return self.data[il] if self.is_cube else self.data


@dataclass
class DiskFitResult:
    model: str
    scale: float
    chisq_binned: float
    nbins: int
    params: dict[str, float]
    psfsub_col: Path
    input_col: Path
    psfsub_cube: Path | None = None
    relatten_cube: Path | None = None
    scale_per_channel: np.ndarray | None = None
    chisq_per_channel: np.ndarray | None = None


def load_disk_model(path: str | Path) -> DiskModel:
    """Model image (first HDU with data) and its GRaTeR header parameters."""
    from astropy.io import fits

    path = Path(path)
    with fits.open(path, memmap=False) as hdul:
        hdu = next((h for h in hdul if h.data is not None), None)
        if hdu is None:
            raise ValueError(f"{path} has no image data")
        data = np.asarray(hdu.data, dtype=float)
        hdr = hdul[0].header.copy()
        hdr.update({k: v for k, v in hdu.header.items() if k in MODEL_KEYS})
    params = {k: float(hdr[k]) for k in MODEL_KEYS if k in hdr}
    return DiskModel(path, data, params)


def roi_mask(
    shape,
    *,
    center=None,
    lrmin: float = 0.0,
    lrmax: float | None = None,
    roirange: tuple[float, float, float] | None = None,
    roioffset: tuple[float, float] = (0.0, 0.0),
    snrmap: np.ndarray | None = None,
    snrlim: float | None = None,
) -> np.ndarray:
    """Boolean evaluation region: radial band ∩ rotated rectangle ∩ SNR > ``snrlim``."""
    center = image_center(shape) if center is None else center
    r = radial_grid(shape, center)
    m = r >= lrmin
    if lrmax is not None:
        m &= r <= lrmax
    if roirange is not None:
        dx, dy, ang = roirange
        m &= roi_rect_rotated(shape, dx, dy, ang, offset=roioffset, center=center)
    if snrmap is not None and snrlim is not None:
        m &= np.isfinite(snrmap) & (snrmap > snrlim)
    if not m.any():
        raise ValueError("the region of interest is empty")
    return m


def scale_model_ls(
    ref: np.ndarray, model: np.ndarray, roi: np.ndarray, sigma: np.ndarray | None = None
) -> float:
    """Least-squares scale factor α minimising Σ_roi ((ref − α·model)/σ)²; ``sigma=None`` is unweighted."""
    w = np.ones_like(ref) if sigma is None else 1.0 / np.where(sigma > 0, sigma, np.nan) ** 2
    sel = roi & np.isfinite(ref) & np.isfinite(model) & np.isfinite(w)
    den = np.sum(w[sel] * model[sel] ** 2)
    if den <= 0:
        raise ValueError("model is zero inside the region of interest")
    return float(np.sum(w[sel] * ref[sel] * model[sel]) / den)


def chisq_map(ref: np.ndarray, model: np.ndarray, sigma: np.ndarray, roi: np.ndarray) -> np.ndarray:
    """Per-pixel χ² inside the ROI (0 elsewhere); pixels with no noise estimate are excluded."""
    out = np.zeros_like(ref)
    sel = roi & np.isfinite(ref) & np.isfinite(model) & np.isfinite(sigma) & (sigma > 0)
    out[sel] = ((ref[sel] - model[sel]) / sigma[sel]) ** 2
    return out


def convolve_model(model: np.ndarray, psf: np.ndarray) -> np.ndarray:
    """Convolve a model image with a (centred, odd-sized) PSF stamp, NaN → 0."""
    m = np.where(np.isfinite(model), model, 0.0)
    return fftconvolve(m, psf / max(psf.sum(), 1e-30) * psf.sum(), mode="same")


def model_in_frame(
    model: np.ndarray,
    rot_deg: float,
    center,
    *,
    h0: Header | None = None,
    dec_deg: float | None = None,
    lat_deg: float | None = None,
    nsc: int = 1,
) -> np.ndarray:
    """North-up model → detector orientation of a frame that needs ``rot_deg`` to be north-up.

    With ``h0``/``dec_deg``/``lat_deg`` and ``nsc > 1`` the model is smeared over
    the exposure (field rotation during the integration), as the planet
    injection does; the IDL code added the same rotated image ``nsc`` times (#186).
    """
    if h0 is not None and dec_deg is not None and lat_deg is not None and nsc > 1:
        dpa = smear_offsets(h0, dec_deg, lat_deg, nsc)
    else:
        dpa = np.zeros(1)
    out = np.zeros_like(model, dtype=float)
    for d in dpa:
        r = rotate_image(model, -(rot_deg + d), center, missing=np.nan)
        out += np.where(np.isfinite(r), r, 0.0) / dpa.size
    return out


def fit_outcomes_line(name: str, params: dict[str, float], scale: float, chisq: float) -> str:
    vals = " ".join(f"{params.get(k, np.nan):10.4f}" for k in MODEL_KEYS)
    return f"{name:30s} {scale:12.5e} {vals} {chisq:10.4f}"


def fit_outcomes_header() -> str:
    return (
        f"{'model':30s} {'scale':>12s} " + " ".join(f"{k:>10s}" for k in MODEL_KEYS) + " chisq_red"
    )
