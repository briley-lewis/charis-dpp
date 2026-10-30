"""Reference-frame selection for ADI / SDI / A-LOCI."""

from __future__ import annotations

import numpy as np

__all__ = ["adi_references", "sdi_references", "most_correlated"]


def adi_references(
    pa: np.ndarray, dtpose: np.ndarray, nf: int, ri: float, fwhm: float, nfwhm: float
) -> np.ndarray:
    """Indices of frames whose field rotation moves a source at ``ri`` by more than ``nfwhm·fwhm`` px.

    The exclusion angle is ``nfwhm·fwhm/ri`` (radians → degrees) plus the smear
    during the target exposure ``dtpose[nf]``.
    """
    dpa = np.abs(pa - pa[nf])
    thresh = np.degrees(nfwhm * fwhm / ri) + dtpose[nf]
    return np.flatnonzero(dpa > thresh)


def sdi_references(wav: np.ndarray, il: int, ri: float, fwhm_il: float, nfwhm: float) -> np.ndarray:
    """Channels whose speckles, rescaled to channel ``il``, move by more than ``nfwhm·fwhm`` px at ``ri``."""
    dr = np.abs(wav[il] / wav - 1.0) * ri
    return np.flatnonzero(dr > nfwhm * fwhm_il)


def most_correlated(cov: np.ndarray, nf: int, indim: np.ndarray, nref: int) -> np.ndarray:
    """Keep the ``nref`` reference frames most correlated with frame ``nf`` (A-LOCI frame selection)."""
    if nref <= 0 or nref >= indim.size:
        return indim
    order = np.argsort(cov[indim, nf])[::-1]
    return indim[order[:nref]]
