"""Readers for the spectral libraries and text tables the pipeline consumes.

Layouts (see ``idl/models``):

* Cruz+2018 templates and Montreal library: FITS primary image ``(3, N)`` =
  wavelength [µm], flux, error.
* Bonnefoy SINFONI library: FITS primary image ``(N, 2)`` = wavelength, flux.
* Vega (``alpha_lyr_stis_005.fits``) and Castelli–Kurucz (``ckp00_*.fits``):
  BINTABLE with ``WAVELENGTH`` [Å] and ``FLUX`` / ``gNN`` columns.
* Pickles ``pickles_uk_N.fits``: BINTABLE ``WAVELENGTH`` [Å], ``FLUX``.
* Two/three-column ASCII spectra (``#`` comments).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

__all__ = [
    "Spectrum",
    "read_spectrum",
    "read_vega",
    "read_kurucz",
    "read_pickles_fits",
    "read_ascii_columns",
]


@dataclass
class Spectrum:
    wavelength_um: np.ndarray
    flux: np.ndarray
    error: np.ndarray | None = None
    units: str = ""  # e.g. "flam", "mJy"
    source: str = ""

    def clip(self, lo_um: float, hi_um: float) -> Spectrum:
        m = (self.wavelength_um >= lo_um) & (self.wavelength_um <= hi_um)
        return Spectrum(
            self.wavelength_um[m],
            self.flux[m],
            None if self.error is None else self.error[m],
            self.units,
            self.source,
        )

    def sorted(self) -> Spectrum:
        o = np.argsort(self.wavelength_um)
        return Spectrum(
            self.wavelength_um[o],
            self.flux[o],
            None if self.error is None else self.error[o],
            self.units,
            self.source,
        )


def read_ascii_columns(path: str | Path, ncols: int | None = None) -> np.ndarray:
    """Numeric columns of a whitespace/comma separated text file; ``#``/``;`` lines are comments."""
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s[0] in "#;":
            continue
        parts = s.replace(",", " ").split()
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        if ncols is not None:
            if len(vals) < ncols:
                continue
            vals = vals[:ncols]
        rows.append(vals)
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows], dtype=float)


def read_spectrum(path: str | Path, *, wavelength_unit: str = "um") -> Spectrum:
    """Read a model/template spectrum in any of the supported layouts; wavelengths returned in µm."""
    path = Path(path)
    if path.suffix.lower() in (".fits", ".fit", ".fts") or path.name.endswith(".fits.gz"):
        with fits.open(path) as hdul:
            data = hdul[0].data
            if data is None:
                return _from_table(Table.read(hdul[1]), str(path))
            data = np.asarray(data, dtype=float)
        if data.ndim != 2:
            raise ValueError(f"{path}: unexpected spectrum shape {data.shape}")
        if data.shape[0] in (2, 3) and data.shape[1] > 3:
            wav, flux = data[0], data[1]
            err = data[2] if data.shape[0] == 3 else None
        elif data.shape[1] in (2, 3):
            wav, flux = data[:, 0], data[:, 1]
            err = data[:, 2] if data.shape[1] == 3 else None
        else:
            raise ValueError(f"{path}: unexpected spectrum shape {data.shape}")
    else:
        arr = read_ascii_columns(path)
        wav, flux = arr[:, 0], arr[:, 1]
        err = arr[:, 2] if arr.shape[1] >= 3 else None
    wav = _to_um(wav, wavelength_unit)
    return Spectrum(wav, flux, err, source=str(path)).sorted()


def _to_um(wav: np.ndarray, unit: str) -> np.ndarray:
    unit = unit.lower()
    if unit in ("um", "micron", "microns"):
        return wav
    if unit in ("nm",):
        return wav * 1e-3
    if unit in ("a", "angstrom", "angstroms"):
        return wav * 1e-4
    raise ValueError(f"unknown wavelength unit {unit!r}")


def _from_table(t: Table, source: str) -> Spectrum:
    cols = {c.upper(): c for c in t.colnames}
    wav = np.asarray(t[cols["WAVELENGTH"]], dtype=float) * 1e-4  # Å → µm
    fcol = cols.get("FLUX", list(cols.values())[1])
    return Spectrum(wav, np.asarray(t[fcol], dtype=float), units="flam", source=source).sorted()


def read_vega(path: str | Path) -> Spectrum:
    """Vega STIS spectrum: wavelength Å → µm, flux erg/s/cm²/Å."""
    return _from_table(Table.read(path, hdu=1), str(path))


def read_kurucz(path: str | Path, logg: float) -> Spectrum:
    """Castelli–Kurucz ``ckp00_<Teff>.fits``: column ``gNN`` nearest to ``logg`` (NN = 10·logg)."""
    t = Table.read(path, hdu=1)
    gcols = [c for c in t.colnames if c.lower().startswith("g") and c[1:].isdigit()]
    if not gcols:
        raise ValueError(f"{path}: no gNN columns")
    want = int(round(logg * 10))
    best = min(gcols, key=lambda c: abs(int(c[1:]) - want))
    wav = np.asarray(t["WAVELENGTH"], dtype=float) * 1e-4
    return Spectrum(
        wav, np.asarray(t[best], dtype=float), units="flam", source=f"{path}[{best}]"
    ).sorted()


def read_pickles_fits(path: str | Path) -> Spectrum:
    return _from_table(Table.read(path, hdu=1), str(path))
