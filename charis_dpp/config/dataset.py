"""Per-dataset configuration (``<target>_<band>.toml``).

Sections (all optional except ``observation.band``/``date`` for resolution):

``[target]``      name, ra_deg, dec_deg, sptype, plx_mas, distance_pc, mags, mag_errors, pm
``[observation]`` date, mjd, band, astrogrid_lambda_d, modulation_nm, frames{science,sky,ladder}
``[geometry]``    frame_shape, star_center, reference_channel, telluric_channels, sub_band_channels
``[reduction]``   the algorithm parameters (today's ``.info`` keys + code-only defaults)
``[pdi]``         HWP log, cycle definition, radii
``[instrument]``  overrides of ``instrument.toml`` keys for this dataset only

:meth:`DatasetConfig.from_info` / :meth:`to_info` convert to and from the
legacy ``.info`` file so the IDL tree can share a dataset directory.
"""

from __future__ import annotations

import copy
import logging
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import tomli_w

from charis_dpp.config.instrument import Instrument, InstrumentTable, load_instrument, to_mjd
from charis_dpp.io.params import ParamFile

log = logging.getLogger(__name__)

SCHEMA_VERSION = 1
TOP_LEVEL = {
    "schema_version",
    "target",
    "observation",
    "geometry",
    "reduction",
    "pdi",
    "instrument",
}

# Defaults for every reduction parameter.  These are the values charis_newobs
# wrote into the .info file plus the keyword defaults that only lived in code.
REDUCTION_DEFAULTS: dict[str, Any] = {
    "fwhm_px": 2.63,
    "rsat_px": -1.0,
    "nfwhm": 0.75,
    "dr_px": 5.0,
    "rmin_px": 5,
    "rmax_px": 65,
    "meanadd": False,
    "zero": False,
    "rsub": False,
    "planet_model": "L0_g.fits",
    "loci": {"na": 100.0, "geom": 1.0, "svd_log10": -6.0, "nref": 50, "pixmask": False},
    "klip": {"npca": 5, "nzone": 1, "meansub": True},
    "register": {
        "method": "sats",
        "smoothval": 3,
        "fittrim": 4,
        "fwhmlim": 7.0,
        "medbox": 15,
        "rin_corr_arcsec": 0.2,
        "rout_corr_arcsec": 0.5,
    },
    "skysub": {"skyrad_px": 60, "skylim": 25.0},
    "specphot": {
        "calmethod": 0,
        "starlib": "pickles",
        "ap_factor": 1.0,
        "av": 0.0,
        "siglim": 5.0,
        "sky_annulus": [2.0, 6.0],
        "sky_annulus_unsat": [4.0, 10.0],
    },
    "covar": {"rho_px": 10.0, "width": 0.5, "step": 0.125, "excval": 3.0},
    "psf": {"fsize_px": 21, "dpsf_px": 21},
    "fwdmod": {
        "dlrstep_px": 2.5,
        "lrmin_px": 10.0,
        "lrmax_px": 60.0,
        "ntc": 5,
        "contrast_log10": -6.0,
        "nsc": 21,
    },
    "detection": {"threshold": 5.0, "excl_px": 2.0, "rangeconstant": 1.5},
}

PDI_DEFAULTS: dict[str, Any] = {
    "hwp_log": "",
    "hwp_source": "log",  # "log" (VAMPIRES csv) or "header" (RET-ANG1)
    "hwp_time_ref": "start",
    "polcent_offset_px": "instrument",
    "astrogrid_lambda_d": "instrument",
    "rlims": {"ip": [10, 25], "bg": [55, 60], "unres": [0, 6]},
    "xcorr": {"upsample": 50},
}


def _deep_merge(base: Mapping, over: Mapping) -> dict:
    out = copy.deepcopy(dict(base))
    for k, v in over.items():
        if isinstance(v, Mapping) and isinstance(out.get(k), Mapping):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _get(d: Mapping, dotted: str, default: Any = None) -> Any:
    cur: Any = d
    for part in dotted.split("."):
        if not isinstance(cur, Mapping) or part not in cur:
            return default
        cur = cur[part]
    return cur


def _set(d: dict, dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    cur = d
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = value


@dataclass
class Resolved:
    """Everything a pipeline step needs, derived once from a :class:`DatasetConfig`."""

    config: DatasetConfig
    instrument: Instrument
    band: str
    wavelengths_um: np.ndarray
    reference_channel: int
    telluric_channels: np.ndarray
    sub_band_channels: dict[str, np.ndarray]
    fwhm_px: np.ndarray  # λ/D per channel, pixels
    frame_shape: tuple[int, int] | None
    star_center: tuple[float, float] | None

    @property
    def nlam(self) -> int:
        return int(self.wavelengths_um.size)

    @property
    def mjd(self) -> float:
        return self.instrument.mjd

    def good_channels(self) -> np.ndarray:
        mask = np.ones(self.nlam, dtype=bool)
        mask[self.telluric_channels] = False
        return np.flatnonzero(mask)


@dataclass
class DatasetConfig:
    data: dict[str, Any] = field(default_factory=dict)
    path: Path | None = None

    # --- construction --------------------------------------------------------------
    @classmethod
    def new(
        cls,
        *,
        target: Mapping[str, Any],
        observation: Mapping[str, Any],
        reduction: Mapping[str, Any] | None = None,
        pdi: Mapping[str, Any] | None = None,
        path: Path | None = None,
    ) -> DatasetConfig:
        data = {
            "schema_version": SCHEMA_VERSION,
            "target": dict(target),
            "observation": dict(observation),
            "geometry": {},
            "reduction": _deep_merge(REDUCTION_DEFAULTS, reduction or {}),
            "pdi": _deep_merge(PDI_DEFAULTS, pdi or {}),
            "instrument": {},
        }
        return cls(data, path)

    @classmethod
    def load(cls, path: str | Path) -> DatasetConfig:
        """Load a TOML dataset file, or a legacy ``.info`` file (converted on the fly)."""
        path = Path(path)
        if path.suffix.lower() == ".info":
            return cls.from_info(path)
        with open(path, "rb") as fh:
            raw = tomllib.load(fh)
        unknown = set(raw) - TOP_LEVEL
        if unknown:
            raise ValueError(f"{path}: unknown top-level sections {sorted(unknown)}")
        raw.setdefault("schema_version", SCHEMA_VERSION)
        raw["reduction"] = _deep_merge(REDUCTION_DEFAULTS, raw.get("reduction", {}))
        raw["pdi"] = _deep_merge(PDI_DEFAULTS, raw.get("pdi", {}))
        for sec in ("target", "observation", "geometry", "instrument"):
            raw.setdefault(sec, {})
        cfg = cls(raw, path)
        cfg.validate()
        return cfg

    def save(self, path: str | Path | None = None, *, write_info: bool = True) -> Path:
        path = Path(path or self.path)
        if path is None:
            raise ValueError("no path given")
        self.validate()
        clean = _strip_none(self.data)
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "wb") as fh:
            tomli_w.dump(clean, fh)
        tmp.replace(path)
        self.path = path
        if write_info:
            self.to_info(path.with_suffix(".info"))
        return path

    # --- access --------------------------------------------------------------------
    def get(self, dotted: str, default: Any = None) -> Any:
        return _get(self.data, dotted, default)

    def set(self, dotted: str, value: Any) -> None:
        _set(self.data, dotted, value)

    def __getitem__(self, section: str) -> dict[str, Any]:
        return self.data[section]

    @property
    def target(self) -> dict[str, Any]:
        return self.data["target"]

    @property
    def observation(self) -> dict[str, Any]:
        return self.data["observation"]

    @property
    def geometry(self) -> dict[str, Any]:
        return self.data["geometry"]

    @property
    def reduction(self) -> dict[str, Any]:
        return self.data["reduction"]

    @property
    def pdi(self) -> dict[str, Any]:
        return self.data["pdi"]

    @property
    def band(self) -> str:
        return str(self.observation.get("band", "lowres"))

    @property
    def mjd(self) -> float:
        obs = self.observation
        if obs.get("mjd") is not None:
            return float(obs["mjd"])
        if obs.get("date"):
            return to_mjd(str(obs["date"]))
        raise ValueError("dataset config needs observation.date or observation.mjd")

    def frames(self, kind: str = "science") -> list[int]:
        from charis_dpp.io.filelists import parse_range

        spec = _get(self.observation, f"frames.{kind}", "")
        return parse_range(spec if isinstance(spec, str) else list(spec))

    def set_frames(self, kind: str, numbers) -> None:
        from charis_dpp.io.filelists import format_range

        _set(self.data, f"observation.frames.{kind}", format_range(numbers))

    @property
    def workdir(self) -> Path:
        return self.path.parent if self.path else Path(".")

    # --- validation -----------------------------------------------------------------
    def validate(self) -> None:
        obs = self.data.get("observation", {})
        band = obs.get("band")
        if band is not None:
            canonical = canonical_band(band)
            if canonical is None:
                raise ValueError(f"observation.band={band!r} is not a known CHARIS band")
            obs["band"] = canonical
        for key in ("date", "mjd"):
            if obs.get(key) is not None:
                to_mjd(obs[key])  # raises on garbage
        shape = self.data.get("geometry", {}).get("frame_shape")
        if shape is not None and (len(shape) != 2 or min(shape) < 8):
            raise ValueError(f"geometry.frame_shape={shape!r} is not a (ny, nx) pair")

    # --- resolution -----------------------------------------------------------------
    def instrument(self, table: InstrumentTable | None = None) -> Instrument:
        table = table or load_instrument()
        return table.for_mjd(self.mjd, self.data.get("instrument") or None)

    def resolve(self, table: InstrumentTable | None = None) -> Resolved:
        from charis_dpp.io.wavelengths import (
            band_channels,
            reference_channel,
            telluric_channels,
            wavelengths_um,
        )

        inst = self.instrument(table)
        band = self.band
        wav = wavelengths_um(band, inst)
        geom = self.geometry
        ref = geom.get("reference_channel", "auto")
        ref = reference_channel(wav, inst, band) if ref in (None, "auto") else int(ref)
        tel = geom.get("telluric_channels", "auto")
        tel = telluric_channels(wav, inst, band) if tel in (None, "auto") else np.asarray(tel, int)
        sub = geom.get("sub_band_channels", "auto")
        sub = (
            band_channels(wav, inst, band)
            if sub in (None, "auto")
            else {k: np.asarray(v, int) for k, v in sub.items()}
        )
        shape = geom.get("frame_shape")
        center = geom.get("star_center")
        return Resolved(
            config=self,
            instrument=inst,
            band=band,
            wavelengths_um=wav,
            reference_channel=ref,
            telluric_channels=tel,
            sub_band_channels=sub,
            fwhm_px=np.asarray(inst.lambda_over_d_px(wav)),
            frame_shape=tuple(int(x) for x in shape) if shape else None,
            star_center=tuple(float(x) for x in center) if center else None,
        )

    # --- legacy .info ---------------------------------------------------------------
    _INFO_TARGET = {
        "NAME": "name",
        "RA": "ra_deg",
        "DEC": "dec_deg",
        "RA_HMS": "ra_hms",
        "DEC_DMS": "dec_dms",
        "PLX": "plx_mas",
        "DISTANCE": "distance_pc",
        "EPLX": "plx_err_mas",
        "PMRA": "pmra_mas_yr",
        "EPMRA": "pmra_err_mas_yr",
        "PMDEC": "pmdec_mas_yr",
        "EPMDEC": "pmdec_err_mas_yr",
        "SPTYPE": "sptype",
    }
    _INFO_MAGS = {"J": ("JMAG", "EJMAG"), "H": ("HMAG", "EHMAG"), "K": ("KMAG", "EKMAG")}
    _INFO_RED = {
        "RSAT": "rsat_px",
        "FWHM": "fwhm_px",
        "NFWHM": "nfwhm",
        "DR": "dr_px",
        "RMIN": "rmin_px",
        "RMAX": "rmax_px",
        "MEANADD": "meanadd",
        "ZERO": "zero",
        "RSUB": "rsub",
        "LOCINA": "loci.na",
        "LOCIGEOM": "loci.geom",
        "ALOCISVD": "loci.svd_log10",
        "ALOCINRE": "loci.nref",
        "PIXMASK": "loci.pixmask",
        "KLIPNPCA": "klip.npca",
        "KLIPNZON": "klip.nzone",
        "MEANSUB": "klip.meansub",
        "PLANMOD": "planet_model",
        "NA": "klip.na",
    }
    _INFO_BOOL = {"MEANADD", "ZERO", "RSUB", "PIXMASK", "MEANSUB"}

    @classmethod
    def from_info(cls, path: str | Path, band: str | None = None) -> DatasetConfig:
        path = Path(path)
        pf = ParamFile.load(path)
        target: dict[str, Any] = {}
        for k, name in cls._INFO_TARGET.items():
            if k in pf and pf[k] not in (None, ""):
                target[name] = pf[k]
        mags, errs = {}, {}
        for b, (mk, ek) in cls._INFO_MAGS.items():
            if mk in pf:
                mags[b] = float(pf[mk])
            if ek in pf:
                errs[b] = float(pf[ek])
        if mags:
            target["mags"] = mags
        if errs:
            target["mag_errors"] = errs
        if band is None:
            stem = path.stem
            band = stem.rsplit("_", 1)[-1] if "_" in stem else "lowres"
        obs: dict[str, Any] = {"band": canonical_band(band) or "lowres"}
        if "OBSDATE" in pf and pf["OBSDATE"]:
            obs["date"] = str(pf["OBSDATE"])
        frames = {}
        for k, kind in (("FNUM_SAT", "science"), ("FNUM_SKY", "sky"), ("FNUM_LAD", "ladder")):
            if k in pf:
                frames[kind] = str(pf[k] or "")
        obs["frames"] = frames
        red: dict[str, Any] = {}
        for k, dotted in cls._INFO_RED.items():
            if k in pf:
                v = pf[k]
                if k in cls._INFO_BOOL:
                    v = bool(int(v))
                _set(red, dotted, v)
        cfg = cls.new(target=target, observation=obs, reduction=red, path=path.with_suffix(".toml"))
        cfg.validate()
        return cfg

    def to_info(self, path: str | Path) -> Path:
        """Write the legacy ``.info`` card file used by the IDL pipeline."""
        pf = ParamFile(path=Path(path))
        t = self.target
        pf.set("name", t.get("name", ""), "")
        for k, name in self._INFO_TARGET.items():
            if k == "NAME":
                continue
            if name in t:
                pf.set(k, t[name], _INFO_COMMENTS.get(k, ""))
        for b, (mk, ek) in self._INFO_MAGS.items():
            if b in t.get("mags", {}):
                pf.set(mk, t["mags"][b], f"{b} magnitude, 2MASS")
            if b in t.get("mag_errors", {}):
                pf.set(ek, t["mag_errors"][b], f"Error on {mk}")
        obs = self.observation
        pf.set("obsdate", str(obs.get("date", "")).replace("-", ""), "Observation date")
        fr = obs.get("frames", {})
        pf.set("fnum_sat", fr.get("science", ""), "File #s, sat, main seq")
        pf.set("fnum_sky", fr.get("sky", ""), "File #s, sky frames")
        pf.set("fnum_lad", fr.get("ladder", ""), "File #s, ladder")
        pf.set("procdir", "./reduc/proc/", "Subdirectory to put processed files")
        r = self.reduction
        pf.set("rsat", float(r["rsat_px"]), "PSF saturation radius")
        pf.set("fwhm", float(r["fwhm_px"]), "Broadband PSF FWHM in pixels")
        pf.set_separator("PSF Sub Range")
        pf.set("nfwhm", float(r["nfwhm"]), "Min displacement for sub, in # of PSF FWHM")
        pf.set("dr", float(r["dr_px"]), "Width of subtraction annuli")
        pf.set("rmin", int(r["rmin_px"]), "Min radius for PSF sub")
        pf.set("rmax", int(r["rmax_px"]), "Max radius for PSF sub")
        pf.set_separator("Switches for Sub/Comb")
        pf.set("meanadd", int(bool(r["meanadd"])), "median-comb(=0), mean-comb(=1)")
        pf.set("zero", int(bool(r["zero"])), "no sub median annulus(=0), sub med. annulus(=1)")
        pf.set("rsub", int(bool(r["rsub"])), "no rad-prof sub.(=0), rad-prof sub. b4 PSF sub(=1)")
        pf.set_separator("(A-)LOCI Params")
        pf.set("locina", float(r["loci"]["na"]), "Area of optimization regions, LOCI")
        pf.set("locigeom", float(r["loci"]["geom"]), "Geometry of optimization regions, LOCI")
        pf.set("alocisvd", r["loci"]["svd_log10"], "log of the SVD cutoff, A-LOCI")
        pf.set("alocinref", int(r["loci"]["nref"]), "frame selection for A-LOCI")
        pf.set("pixmask", int(bool(r["loci"]["pixmask"])), "no pixmask(=0),pixmask(=1)")
        pf.set_separator("KLIP Params")
        pf.set("klipnpca", int(r["klip"]["npca"]), "Number of KL Modes, KLIP")
        pf.set("klipnzone", int(r["klip"]["nzone"]), "Number of azimuthal zones per annulus, KLIP")
        pf.set(
            "meansub", int(bool(r["klip"]["meansub"])), "no subtract mean(=0), subtract mean(=1)"
        )
        pf.set_separator("PlntFwdMod")
        pf.set("planmod", str(r["planet_model"]), "Model Plan Spec")
        return pf.save()


_INFO_COMMENTS = {
    "RA": "RA, equinox J2000, epoch J2000",
    "DEC": "DEC, equinox J2000, epoch J2000",
    "RA_HMS": "RA, equinox J2000, epoch J2000",
    "DEC_DMS": "DEC, equinox J2000, epoch J2000",
    "PLX": " Parallax, mas",
    "DISTANCE": " Distance (pc)",
    "EPLX": " Error on parallax, mas",
    "PMRA": " Proper motion RA [*cos(DEC)], mas/yr",
    "EPMRA": " Error on PMRA",
    "PMDEC": " Proper motion DEC, mas/yr",
    "EPMDEC": " Error on PMDEC",
    "SPTYPE": " Spectral Type",
}


def _strip_none(d: Any) -> Any:
    if isinstance(d, Mapping):
        return {k: _strip_none(v) for k, v in d.items() if v is not None}
    if isinstance(d, list | tuple):
        return [_strip_none(v) for v in d]
    if isinstance(d, np.generic):
        return d.item()
    if isinstance(d, np.ndarray):
        return d.tolist()
    return d


def canonical_band(name: str | None) -> str | None:
    """Map any CHARIS band spelling (``Broadband``, ``low``, ``j``…) to ``lowres``/``J``/``H``/``K``."""
    if name is None:
        return None
    n = str(name).strip()
    table = load_instrument()
    for band, spec in table.bands.items():
        if n == band or n in spec.get("aliases", []):
            return band
    low = n.lower()
    for band, spec in table.bands.items():
        if low == band.lower() or low in [a.lower() for a in spec.get("aliases", [])]:
            return band
    return None
