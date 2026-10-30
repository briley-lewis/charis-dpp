"""CHARIS cube I/O.

Pipeline products are two-HDU FITS files: an empty primary HDU carrying the
observation header (``h0``) and an image extension holding the cube with its
own header (``h1``).  Data are returned as ``(nλ, ny, nx)`` float arrays —
astropy's natural order for an IDL ``[x, y, λ]`` cube — so code indexes
``cube[l, y, x]``.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.io.fits import Header


@dataclass
class Cube:
    data: np.ndarray
    h0: Header = field(default_factory=Header)
    h1: Header = field(default_factory=Header)
    path: Path | None = None

    @property
    def nlam(self) -> int:
        return int(self.data.shape[0]) if self.data.ndim == 3 else 1

    @property
    def shape(self) -> tuple[int, int]:
        return int(self.data.shape[-2]), int(self.data.shape[-1])

    def copy(self) -> Cube:
        return Cube(self.data.copy(), self.h0.copy(), self.h1.copy(), self.path)

    def __getitem__(self, key: str) -> Any:
        """Header lookup: extension header first, then primary."""
        if key in self.h1:
            return self.h1[key]
        return self.h0[key]

    def get(self, key: str, default: Any = None) -> Any:
        if key in self.h1:
            return self.h1[key]
        return self.h0.get(key, default)

    def set(self, key: str, value: Any, comment: str | None = None, both: bool = False) -> None:
        self.h1[key] = (value, comment) if comment else value
        if both:
            self.h0[key] = (value, comment) if comment else value


def read_cube(path: str | Path, ext: int = 1, dtype=np.float64) -> Cube:
    """Read a pipeline cube (or a single-HDU file) as a :class:`Cube`."""
    path = Path(path)
    with fits.open(path, memmap=False) as hdul:
        h0 = hdul[0].header.copy()
        if len(hdul) > ext and hdul[ext].data is not None:
            data = np.asarray(hdul[ext].data, dtype=dtype)
            h1 = hdul[ext].header.copy()
        elif hdul[0].data is not None:
            data = np.asarray(hdul[0].data, dtype=dtype)
            h1 = h0.copy()
        else:
            raise ValueError(f"{path}: no image data in HDU 0 or {ext}")
    return Cube(data, h0, h1, path)


def read_header(path: str | Path, ext: int = 0) -> Header:
    return fits.getheader(path, ext=ext)


def read_image(path: str | Path, ext: int | None = None, dtype=np.float64) -> np.ndarray:
    """Image/cube data only; takes the first HDU with data unless ``ext`` is given."""
    with fits.open(path, memmap=False) as hdul:
        if ext is not None:
            return np.asarray(hdul[ext].data, dtype=dtype)
        for hdu in hdul:
            if hdu.data is not None:
                return np.asarray(hdu.data, dtype=dtype)
    raise ValueError(f"{path}: no image data")


def write_cube(
    path: str | Path,
    data: np.ndarray | Cube,
    h0: Header | None = None,
    h1: Header | None = None,
    *,
    dtype=np.float32,
    overwrite: bool = True,
) -> Path:
    """Write the standard two-HDU layout (empty primary + image extension)."""
    path = Path(path)
    if isinstance(data, Cube):
        h0 = data.h0 if h0 is None else h0
        h1 = data.h1 if h1 is None else h1
        data = data.data
    h0 = Header() if h0 is None else h0.copy()
    h1 = Header() if h1 is None else h1.copy()
    for h in (h0, h1):
        for k in (
            "NAXIS1",
            "NAXIS2",
            "NAXIS3",
            "BITPIX",
            "SIMPLE",
            "EXTEND",
            "XTENSION",
            "PCOUNT",
            "GCOUNT",
        ):
            h.pop(k, None)
    primary = fits.PrimaryHDU(header=h0)
    ext = fits.ImageHDU(data=np.asarray(data, dtype=dtype), header=h1)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([primary, ext]).writeto(path, overwrite=overwrite)
    return path


def write_image(
    path: str | Path, data: np.ndarray, header: Header | None = None, *, dtype=np.float32
) -> Path:
    """Single-HDU file (attenuation maps, PSF cubes, diagnostics)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=np.asarray(data, dtype=dtype), header=header).writeto(path, overwrite=True)
    return path


def header_get(header: Header, key: str, default: Any = None, *, justfirst: bool = True) -> Any:
    """``sxpar`` with ``/justfirst``: the first of duplicated keywords (RA/DEC appear twice)."""
    if key not in header:
        return default
    values = header[key] if not isinstance(header[key], list | tuple) else header[key]
    try:
        all_vals = [c.value for c in header.cards if c.keyword == key.upper()]
    except Exception:  # pragma: no cover
        return values
    if not all_vals:
        return default
    return all_vals[0] if justfirst else all_vals[-1]


# --- satellite spots --------------------------------------------------------------
def get_sat_positions(h1: Header, nlam: int | None = None) -> np.ndarray:
    """Satellite-spot centroids from ``SATS<l>_<j>`` cards → ``(nλ, 4, 2)`` of (x, y); NaN if missing."""
    if nlam is None:
        nlam = int(h1.get("NAXIS3", 0)) or _count_sats(h1)
    out = np.full((nlam, 4, 2), np.nan)
    for l in range(nlam):
        for j in range(4):
            key = f"SATS{l}_{j}"
            if key in h1:
                parts = str(h1[key]).split()
                if len(parts) >= 2:
                    out[l, j] = float(parts[0]), float(parts[1])
    return out


def _count_sats(h1: Header) -> int:
    n = 0
    while f"SATS{n}_0" in h1:
        n += 1
    return n


def set_sat_positions(h1: Header, sats: np.ndarray) -> None:
    for l in range(sats.shape[0]):
        for j in range(sats.shape[1]):
            x, y = sats[l, j]
            h1[f"SATS{l}_{j}"] = (f"{x:.3f} {y:.3f}", "sat spot centroid x y (0-based)")


def get_psf_center(
    h1: Header, default: tuple[float, float] | None = None
) -> tuple[float, float] | None:
    """Star position written by registration (``PSFCENTX/Y``)."""
    if "PSFCENTX" in h1 and "PSFCENTY" in h1:
        return float(h1["PSFCENTX"]), float(h1["PSFCENTY"])
    return default


def copy_cards(src: Header, dst: Header, keys: Iterable[str]) -> None:
    for k in keys:
        if k in src:
            dst[k] = (src[k], src.comments[k])
