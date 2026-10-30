"""Instrument constants resolved by observation date.

``instrument.toml`` holds dated ``[[epoch]]`` blocks; :func:`load_instrument`
parses the file and :meth:`InstrumentTable.for_mjd` merges every epoch whose
``valid_from`` is on or before the requested MJD, later epochs overriding
earlier ones key by key.  The result is an immutable :class:`Instrument`.

Replaces ``idl/setup/charis_get_constant.pro`` and the constants scattered
through ``get_charis_sat_contrast.pro``, ``charis_imprep.pro``,
``charis_nbs_corr.pro`` and ``pdi/**``.
"""

from __future__ import annotations

import copy
import logging
import os
import tomllib
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from datetime import date, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from astropy.time import Time

from charis_dpp.config.paths import package_data

log = logging.getLogger(__name__)

_MISSING = object()


class Section(Mapping[str, Any]):
    """Read-only mapping with attribute access and dotted ``get``."""

    __slots__ = ("_d",)

    def __init__(self, d: Mapping[str, Any]):
        object.__setattr__(self, "_d", {k: _wrap(v) for k, v in d.items()})

    def __getattr__(self, name: str) -> Any:
        try:
            return self._d[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __getitem__(self, key: str) -> Any:
        return self._d[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._d)

    def __len__(self) -> int:
        return len(self._d)

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("Section is read-only")

    def get(self, key: str, default: Any = None) -> Any:  # type: ignore[override]
        """``get("a.b.c")`` walks nested sections."""
        cur: Any = self
        for part in key.split("."):
            if not isinstance(cur, Section) or part not in cur._d:
                return default
            cur = cur._d[part]
        return cur

    def to_dict(self) -> dict[str, Any]:
        return _unwrap(self)

    def __repr__(self) -> str:
        return f"Section({self._d!r})"


def _wrap(v: Any) -> Any:
    if isinstance(v, Mapping) and not isinstance(v, Section):
        return Section(v)
    if isinstance(v, list):
        return tuple(_wrap(x) for x in v)
    return v


def _unwrap(v: Any) -> Any:
    if isinstance(v, Section):
        return {k: _unwrap(x) for k, x in v._d.items()}
    if isinstance(v, tuple):
        return [_unwrap(x) for x in v]
    return v


def _deep_merge(base: dict, over: Mapping) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, Mapping) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def to_mjd(when: float | int | str | date | datetime | Time) -> float:
    """Accept an MJD, ISO date/datetime string, ``date``/``datetime`` or ``Time``."""
    if isinstance(when, Time):
        return float(when.mjd)
    if isinstance(when, int | float | np.floating | np.integer):
        return float(when)
    if isinstance(when, datetime):
        return float(Time(when).mjd)
    if isinstance(when, date):
        return float(Time(datetime(when.year, when.month, when.day)).mjd)
    if isinstance(when, str):
        s = when.strip()
        if len(s) == 8 and s.isdigit():  # CHARIS 'YYYYMMDD'
            s = f"{s[:4]}-{s[4:6]}-{s[6:]}"
        return float(Time(s, format="isot" if "T" in s else "iso").mjd)
    raise TypeError(f"cannot convert {when!r} to MJD")


@dataclass(frozen=True)
class Epoch:
    label: str
    valid_from_mjd: float
    source: str
    values: dict[str, Any]
    placeholder_date: bool = False


@dataclass(frozen=True)
class Instrument:
    """Constants in force at one observation date (immutable)."""

    mjd: float
    epoch_label: str
    epochs_applied: tuple[str, ...]
    site: Section
    detector: Section
    bands: Section
    photometry: Section
    values: Section  # merged epoch block
    overrides: dict[str, Any] = field(default_factory=dict)

    # --- frequently used scalars -------------------------------------------------
    @property
    def dtel_m(self) -> float:
        return float(self.values.dtel_m)

    @property
    def pixscale_arcsec(self) -> float:
        return float(self.values.pixscale_arcsec)

    @property
    def north_offset_deg(self) -> float:
        return float(self.values.north_offset_deg)

    @property
    def parity(self) -> int:
        return int(self.values.parity)

    @property
    def post_nbs(self) -> bool:
        return self.parity == -1

    @property
    def latitude_deg(self) -> float:
        return float(self.site.latitude_deg)

    @property
    def longitude_deg(self) -> float:
        return float(self.site.longitude_deg)

    @property
    def read_time_s(self) -> float:
        return float(self.detector.read_time_s)

    @property
    def sat_spots(self) -> Section:
        return self.values.sat_spots

    @property
    def pdi(self) -> Section:
        return self.values.pdi

    def get(self, key: str, default: Any = _MISSING) -> Any:
        """Dotted lookup over the merged epoch values, e.g. ``"pdi.mueller.m3.m1"``."""
        v = self.values.get(key, _MISSING)
        if v is _MISSING:
            if default is _MISSING:
                raise KeyError(f"instrument constant {key!r} not defined for MJD {self.mjd:.1f}")
            return default
        return v

    def lambda_over_d_px(self, wavelength_um: np.ndarray | float) -> np.ndarray | float:
        """λ/D in pixels (the IDL per-channel ``fwhm``)."""
        rad_to_arcsec = 206264.806
        return np.asarray(wavelength_um) * 1e-6 / self.dtel_m * rad_to_arcsec / self.pixscale_arcsec

    def spot_guesses_px(self, band: str, astrogrid_lambda_d: float | None = None) -> np.ndarray:
        """Default satellite-spot (x, y) offsets from the frame centre, (4, 2)."""
        key = band if band in self.sat_spots.guess_xy else "lowres"
        g = np.array(self.sat_spots.guess_xy[key], dtype=float)
        nominal = float(self.sat_spots.nominal_astrogrid_lambda_d)
        if astrogrid_lambda_d is not None:
            g = g * (astrogrid_lambda_d / nominal)
        if self.parity == -1:
            g = g * np.array([-1.0, 1.0])
        return g

    def to_header_cards(self) -> list[tuple[str, Any, str]]:
        """Cards recording the constants used, for output FITS headers."""
        return [
            ("INSTEPOC", self.epoch_label, "instrument epoch (charis_dpp instrument.toml)"),
            ("INSTMJD", round(self.mjd, 5), "MJD used to resolve instrument constants"),
            ("DTEL", self.dtel_m, "telescope diameter [m]"),
            ("PIXSCALE", self.pixscale_arcsec, "pixel scale [arcsec/px]"),
            ("NORTHOFF", self.north_offset_deg, "north offset applied for north-up [deg]"),
            ("PARITY", self.parity, "+1 as delivered, -1 x-flipped (post-NBS)"),
        ]


@dataclass(frozen=True)
class InstrumentTable:
    """Parsed ``instrument.toml``."""

    path: Path
    schema_version: int
    site: dict[str, Any]
    detector: dict[str, Any]
    bands: dict[str, Any]
    photometry: dict[str, Any]
    epochs: tuple[Epoch, ...]

    def for_mjd(self, when: Any, overrides: Mapping[str, Any] | None = None) -> Instrument:
        mjd = to_mjd(when)
        applied: list[str] = []
        merged: dict[str, Any] = {}
        for ep in self.epochs:
            if ep.valid_from_mjd <= mjd:
                merged = _deep_merge(merged, ep.values)
                applied.append(ep.label)
                if ep.placeholder_date and abs(mjd - ep.valid_from_mjd) < 92:
                    log.warning(
                        "instrument epoch %r has a placeholder start date (%s); the observation "
                        "(MJD %.1f) is within 3 months of it — confirm which side it falls on",
                        ep.label,
                        ep.source,
                        mjd,
                    )
        if not applied:
            raise ValueError(
                f"MJD {mjd:.1f} is before the first instrument epoch "
                f"({self.epochs[0].valid_from_mjd:.1f})"
            )
        ov = dict(overrides or {})
        if ov:
            merged = _deep_merge(merged, ov)
            log.info("instrument overrides applied: %s", ov)
        return Instrument(
            mjd=mjd,
            epoch_label=applied[-1],
            epochs_applied=tuple(applied),
            site=Section(self.site),
            detector=Section(self.detector),
            bands=Section(self.bands),
            photometry=Section(self.photometry),
            values=Section(merged),
            overrides=ov,
        )

    def epoch_labels(self) -> list[str]:
        return [e.label for e in self.epochs]


def _parse(path: Path) -> InstrumentTable:
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    epochs: list[Epoch] = []
    for block in raw.get("epoch", []):
        block = dict(block)
        label = block.pop("label")
        source = block.pop("source", "")
        placeholder = bool(block.pop("placeholder_date", False))
        vf = block.pop("valid_from")
        mjd = float(block.pop("valid_from_mjd", to_mjd(vf)))
        epochs.append(Epoch(label, mjd, source, block, placeholder))
    epochs.sort(key=lambda e: e.valid_from_mjd)
    return InstrumentTable(
        path=path,
        schema_version=int(raw.get("schema_version", 1)),
        site=raw["site"],
        detector=raw["detector"],
        bands=raw["bands"],
        photometry=raw["photometry"],
        epochs=tuple(epochs),
    )


@lru_cache(maxsize=4)
def _load_cached(path_str: str) -> InstrumentTable:
    return _parse(Path(path_str))


def load_instrument(path: str | os.PathLike | None = None) -> InstrumentTable:
    """Load ``instrument.toml`` (package default, ``$CHARIS_DPP_INSTRUMENT``, or ``path``)."""
    if path is None:
        env = os.environ.get("CHARIS_DPP_INSTRUMENT")
        path = Path(env).expanduser() if env else package_data("instrument.toml")
    return _load_cached(str(Path(path).resolve()))


def instrument_for(when: Any, overrides: Mapping[str, Any] | None = None) -> Instrument:
    """Convenience: ``load_instrument().for_mjd(when, overrides)``."""
    return load_instrument().for_mjd(when, overrides)
