"""Time and coordinate helpers (``parangle.pro``, ``rot_ratef.pro``, ``ct2lst``, ``precess``, ``ten``/``sixty``).

Angles are in degrees unless stated.  Hour angles are in hours, like the IDL
code and the ``HA`` header keyword.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.coordinates import FK5, SkyCoord
from astropy.time import Time

__all__ = [
    "ten",
    "sixty",
    "parallactic_angle",
    "rotation_rate",
    "lst_hours",
    "hour_angle",
    "precess",
    "wrap180",
]


def ten(value) -> float:
    """Sexagesimal (``'hh:mm:ss'``, ``'dd mm ss'``, or a number) → decimal."""
    if isinstance(value, int | float | np.floating | np.integer):
        return float(value)
    s = str(value).strip().replace(":", " ")
    parts = [float(p) for p in s.split()]
    if not parts:
        raise ValueError(f"cannot parse {value!r}")
    sign = -1.0 if s.lstrip().startswith("-") else 1.0
    parts[0] = abs(parts[0])
    total = 0.0
    for i, p in enumerate(parts[:3]):
        total += p / 60.0**i
    return sign * total


def sixty(value: float) -> tuple[float, float, float]:
    sign = -1 if value < 0 else 1
    v = abs(value)
    d = int(v)
    m = int((v - d) * 60)
    s = (v - d - m / 60) * 3600
    return (sign * d, m, s)


def wrap180(a):
    return (np.asarray(a, dtype=float) + 180.0) % 360.0 - 180.0


def parallactic_angle(ha_hours, dec_deg: float, lat_deg: float):
    """Parallactic angle (deg, in (-180, 180]) from hour angle (hours), dec and latitude."""
    ha = np.deg2rad(np.asarray(ha_hours, dtype=float) * 15.0)
    dec = np.deg2rad(dec_deg)
    lat = np.deg2rad(lat_deg)
    num = np.sin(ha) * np.cos(lat)
    den = np.sin(lat) * np.cos(dec) - np.cos(lat) * np.sin(dec) * np.cos(ha)
    return np.degrees(np.arctan2(num, den))


def rotation_rate(ha_hours, dec_deg: float, lat_deg: float):
    """Field-rotation rate dPA/dt in deg/s (``rot_ratef``, with an unambiguous azimuth)."""
    ha = np.deg2rad(np.asarray(ha_hours, dtype=float) * 15.0)
    dec = np.deg2rad(dec_deg)
    lat = np.deg2rad(lat_deg)
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(ha)
    alt = np.arcsin(np.clip(sin_alt, -1, 1))
    az = np.arctan2(
        -np.sin(ha) * np.cos(dec),
        np.cos(lat) * np.sin(dec) - np.sin(lat) * np.cos(dec) * np.cos(ha),
    )
    omega = 2 * np.pi / 86164.0905  # sidereal rate, rad/s
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = omega * np.cos(lat) * np.cos(az) / np.cos(alt)
    return np.degrees(rate)


def lst_hours(mjd, lng_deg: float) -> float:
    """Local mean sidereal time in hours (astrolib ``ct2lst`` is mean sidereal time)."""
    t = Time(np.asarray(mjd, dtype=float), format="mjd", scale="utc")
    return t.sidereal_time("mean", longitude=lng_deg * u.deg).hour


def hour_angle(lst_h, ra_deg: float) -> float:
    """HA (hours) in [-12, 12)."""
    ha = (np.asarray(lst_h, dtype=float) - ra_deg / 15.0 + 12.0) % 24.0 - 12.0
    return ha


def precess(
    ra_deg: float, dec_deg: float, equinox_from: float, equinox_to: float
) -> tuple[float, float]:
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame=FK5(equinox=f"J{equinox_from}"))
    c2 = c.transform_to(FK5(equinox=f"J{equinox_to}"))
    return float(c2.ra.deg), float(c2.dec.deg)


def jd_from_date_time(date_iso: str, utc_hours: float) -> float:
    t = Time(date_iso, format="iso", scale="utc") + utc_hours * u.hour
    return float(t.jd)
