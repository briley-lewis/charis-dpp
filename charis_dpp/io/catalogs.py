"""Target information from SIMBAD / Hipparcos / 2MASS (``charis_newobs`` queries), with an offline path.

:class:`StarInfo` is what ``newobs`` writes into the dataset config.  Tests
and offline use construct it directly; :func:`query_star` fills it from the
network with astroquery.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

log = logging.getLogger(__name__)


@dataclass
class StarInfo:
    name: str
    ra_deg: float
    dec_deg: float
    sptype: str = ""
    plx_mas: float = np.nan
    plx_err_mas: float = np.nan
    pmra_mas_yr: float = np.nan
    pmra_err_mas_yr: float = np.nan
    pmdec_mas_yr: float = np.nan
    pmdec_err_mas_yr: float = np.nan
    mags: dict[str, float] = field(default_factory=dict)
    mag_errors: dict[str, float] = field(default_factory=dict)

    @property
    def distance_pc(self) -> float:
        return 1e3 / self.plx_mas if np.isfinite(self.plx_mas) and self.plx_mas > 0 else np.nan

    @property
    def ra_hms(self) -> str:
        from astropy import units as u
        from astropy.coordinates import Angle

        return Angle(self.ra_deg * u.deg).to_string(unit=u.hour, sep=":", precision=2, pad=True)

    @property
    def dec_dms(self) -> str:
        from astropy import units as u
        from astropy.coordinates import Angle

        return Angle(self.dec_deg * u.deg).to_string(
            unit=u.deg, sep=":", precision=1, pad=True, alwayssign=True
        )

    def to_config(self) -> dict[str, Any]:
        d = asdict(self)
        d["distance_pc"] = self.distance_pc
        d["ra_hms"] = self.ra_hms
        d["dec_dms"] = self.dec_dms
        return {k: v for k, v in d.items() if not (isinstance(v, float) and np.isnan(v))}

    @classmethod
    def from_config(cls, target: dict[str, Any]) -> StarInfo:
        keys = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in target.items() if k in keys})


def _first(tbl, col: str, default=np.nan):
    if tbl is None or len(tbl) == 0 or col not in tbl.colnames:
        return default
    v = tbl[col][0]
    try:
        if np.ma.is_masked(v):
            return default
        return float(v) if not isinstance(v, str | bytes) else v
    except (TypeError, ValueError):
        return default


def query_star(name: str, radius_arcmin: float = 20.0 / 60.0, timeout: float = 60.0) -> StarInfo:
    """Resolve ``name`` with SIMBAD, then Hipparcos (I/239) and 2MASS (II/246) via VizieR.

    Tries the name variants ``charis_newobs`` tried (``* name``, ``V* name``).
    Missing values stay NaN (the IDL code wrote 0 or -999 inconsistently).
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier

    Simbad.TIMEOUT = timeout
    sim = Simbad()
    try:
        sim.add_votable_fields("sp_type", "plx_value", "plx_err", "pmra", "pmdec")
    except Exception:  # field names vary across astroquery versions
        pass
    result = None
    for variant in (name, name.replace(" ", ""), f"* {name}", f"V* {name}"):
        try:
            result = sim.query_object(variant)
        except Exception as exc:  # pragma: no cover - network
            log.warning("SIMBAD query %r failed: %s", variant, exc)
            result = None
        if result is not None and len(result) > 0:
            break
    if result is None or len(result) == 0:
        raise LookupError(f"SIMBAD cannot resolve {name!r}")
    cols = {c.lower(): c for c in result.colnames}
    ra_col = cols.get("ra", cols.get("ra_d"))
    dec_col = cols.get("dec", cols.get("dec_d"))
    ra_raw, dec_raw = result[ra_col][0], result[dec_col][0]
    if isinstance(ra_raw, str | bytes):
        c = SkyCoord(f"{ra_raw} {dec_raw}", unit=(u.hourangle, u.deg))
    else:
        c = SkyCoord(float(ra_raw) * u.deg, float(dec_raw) * u.deg)
    info = StarInfo(name=name, ra_deg=float(c.ra.deg), dec_deg=float(c.dec.deg))
    sp = _first(result, cols.get("sp_type", "sp_type"), "")
    info.sptype = str(sp) if sp else ""

    viz = Vizier(columns=["**"], row_limit=5)
    viz.TIMEOUT = timeout
    try:
        hip = viz.query_region(c, radius=radius_arcmin * u.arcmin, catalog="I/239/hip_main")
        hip = hip[0] if len(hip) else None
    except Exception as exc:  # pragma: no cover - network
        log.warning("Hipparcos query failed: %s", exc)
        hip = None
    if hip is not None:
        info.plx_mas = _first(hip, "Plx")
        info.plx_err_mas = _first(hip, "e_Plx")
        info.pmra_mas_yr = _first(hip, "pmRA")
        info.pmra_err_mas_yr = _first(hip, "e_pmRA")
        info.pmdec_mas_yr = _first(hip, "pmDE")
        info.pmdec_err_mas_yr = _first(hip, "e_pmDE")
        if not info.sptype:
            info.sptype = str(_first(hip, "SpType", "") or "")
    else:
        # SIMBAD fallback for parallax / proper motion
        info.plx_mas = _first(result, cols.get("plx_value", "plx_value"))
        info.plx_err_mas = _first(result, cols.get("plx_err", "plx_err"))
        info.pmra_mas_yr = _first(result, cols.get("pmra", "pmra"))
        info.pmdec_mas_yr = _first(result, cols.get("pmdec", "pmdec"))
    try:
        tm = viz.query_region(c, radius=radius_arcmin * u.arcmin, catalog="II/246/out")
        tm = tm[0] if len(tm) else None
    except Exception as exc:  # pragma: no cover - network
        log.warning("2MASS query failed: %s", exc)
        tm = None
    if tm is not None and len(tm):
        if "_r" in tm.colnames:
            tm.sort("_r")
        for band, col in (("J", "Jmag"), ("H", "Hmag"), ("K", "Kmag")):
            v = _first(tm, col)
            e = _first(tm, f"e_{col}")
            if np.isfinite(v):
                info.mags[band] = v
            if np.isfinite(e):
                info.mag_errors[band] = e
    else:
        log.warning("no 2MASS entry within %.1f arcmin of %s", radius_arcmin, name)
    return info
