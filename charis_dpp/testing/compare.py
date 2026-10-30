"""Numerical comparison of two FITS products or two ``reduc.log`` files.

Used by the regression tests (``tests/regression``) against archived IDL
outputs and by ``charis compare-fits``.  The interior of each image (a border
of ``border`` pixels is ignored, where the IDL cubic interpolation edge
behaviour differs) is compared with a relative RMS, a Pearson correlation and
the fraction of finite pixels that agree.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.io.fits import read_image
from charis_dpp.io.reduclog import ReducLog

__all__ = ["CompareResult", "compare_fits", "compare_reduclog"]


@dataclass
class CompareResult:
    name: str
    rel_rms: float
    pearson: float
    finite_agree: float
    max_abs: float
    n: int

    def ok(
        self, *, rel_rms: float = 1e-3, pearson: float = 0.995, finite_agree: float = 0.98
    ) -> bool:
        return (
            self.rel_rms <= rel_rms
            and self.pearson >= pearson
            and self.finite_agree >= finite_agree
        )

    def __str__(self) -> str:
        return (
            f"{self.name}: rel RMS {self.rel_rms:.3e}, r = {self.pearson:.5f}, "
            f"finite agreement {100 * self.finite_agree:.1f}%, max |diff| {self.max_abs:.3e} ({self.n} px)"
        )


def compare_arrays(
    a: np.ndarray, b: np.ndarray, *, border: int = 2, name: str = ""
) -> CompareResult:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"{name}: shapes differ {a.shape} vs {b.shape}")
    if border > 0 and a.ndim >= 2:
        sl = (Ellipsis, slice(border, -border), slice(border, -border))
        a, b = a[sl], b[sl]
    fa, fb = np.isfinite(a), np.isfinite(b)
    agree = float(np.mean(fa == fb)) if a.size else 1.0
    both = fa & fb
    if both.sum() < 2:
        return CompareResult(name, np.inf, 0.0, agree, np.inf, int(both.sum()))
    x, y = a[both], b[both]
    diff = x - y
    scale = np.sqrt(np.mean(x**2)) or 1.0
    rel = float(np.sqrt(np.mean(diff**2)) / scale)
    r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 1.0
    return CompareResult(name, rel, r, agree, float(np.max(np.abs(diff))), int(both.sum()))


def compare_fits(
    a: str | Path, b: str | Path, *, border: int = 2, ext: int | None = None
) -> CompareResult:
    """Compare the image data of two FITS files (first HDU with data unless ``ext``)."""
    return compare_arrays(
        read_image(a, ext=ext),
        read_image(b, ext=ext),
        border=border,
        name=f"{Path(a).name} vs {Path(b).name}",
    )


def compare_reduclog(a: str | Path, b: str | Path) -> dict[str, float]:
    """Max absolute differences of the centre, HA and PA columns for the common frames."""
    la, lb = ReducLog.load(a), ReducLog.load(b)
    common = sorted(set(la.entries) & set(lb.entries))
    if not common:
        raise ValueError("no common frame numbers")
    ca, cb = la.centers(common), lb.centers(common)
    return {
        "n": float(len(common)),
        "max_dx": float(np.max(np.abs(ca[:, 0] - cb[:, 0]))),
        "max_dy": float(np.max(np.abs(ca[:, 1] - cb[:, 1]))),
        "max_dha": float(np.max(np.abs(la.ha(common) - lb.ha(common)))),
        "max_dpa": float(np.max(np.abs(la.pa(common) - lb.pa(common)))),
    }
