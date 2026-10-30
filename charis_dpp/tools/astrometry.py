"""Background-star track: where a stationary background object would appear relative
to the star over time (``tools/astrometry/charis_est_prop_mot.pro``).

The star's proper motion and parallax are applied with astropy
(``SkyCoord.apply_space_motion`` plus the parallactic ellipse from the
observer's position), so the offsets are in the equatorial (E, N) frame
(the IDL mixed ecliptic offsets with equatorial proper motion, #242).  A
position file ``e n sig_e sig_n date`` compares measured epochs with the
prediction (#243).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time

from charis_dpp.io.catalogs import StarInfo

log = logging.getLogger(__name__)

__all__ = ["BackgroundTrack", "background_track", "read_position_file"]


@dataclass
class BackgroundTrack:
    times: Time
    d_east_arcsec: np.ndarray  # offset of a background object relative to the star, E positive
    d_north_arcsec: np.ndarray
    predicted_at_epochs: np.ndarray | None = None  # (nepoch, 2)
    chisq: float | None = None


def read_position_file(
    path: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    rows = [
        ln.split()
        for ln in Path(path).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    e = np.array([float(r[0]) for r in rows])
    n = np.array([float(r[1]) for r in rows])
    se = np.array([float(r[2]) for r in rows])
    sn = np.array([float(r[3]) for r in rows])
    dates = [r[4] for r in rows]
    return e, n, se, sn, dates


def _apparent(star: StarInfo, t: Time) -> SkyCoord:
    """ICRS position moved by proper motion plus the annual parallax ellipse (no aberration:
    it is common to the star and a background object and cancels in relative astrometry)."""
    plx = star.plx_mas if np.isfinite(star.plx_mas) and star.plx_mas > 0 else 0.0
    c = SkyCoord(
        ra=star.ra_deg * u.deg,
        dec=star.dec_deg * u.deg,
        distance=(1000.0 / max(plx, 1e-3)) * u.pc,
        pm_ra_cosdec=(star.pmra_mas_yr if np.isfinite(star.pmra_mas_yr) else 0.0) * u.mas / u.yr,
        pm_dec=(star.pmdec_mas_yr if np.isfinite(star.pmdec_mas_yr) else 0.0) * u.mas / u.yr,
        obstime=Time("J2000"),
        frame="icrs",
    )
    moved = c.apply_space_motion(new_obstime=t)
    earth = get_body_barycentric("earth", t)  # AU
    x, y, z = (earth.x.to(u.au).value, earth.y.to(u.au).value, earth.z.to(u.au).value)
    ra, dec = moved.ra.rad, moved.dec.rad
    p_as = plx * 1e-3  # parallax in arcsec
    d_ra_cosdec = p_as * (x * np.sin(ra) - y * np.cos(ra))
    d_dec = p_as * (x * np.cos(ra) * np.sin(dec) + y * np.sin(ra) * np.sin(dec) - z * np.cos(dec))
    return SkyCoord(
        ra=moved.ra + (d_ra_cosdec / np.cos(dec)) * u.arcsec,
        dec=moved.dec + d_dec * u.arcsec,
        frame="icrs",
    )


def background_track(
    star: StarInfo,
    start_date: str,
    *,
    nyears: float | None = None,
    ndays: int | None = None,
    position: tuple[float, float] = (0.0, 0.0),
    posfile: str | Path | None = None,
    plot: str | Path | None = None,
) -> BackgroundTrack:
    """Offsets (E, N) in arcsec of a fixed background object starting at ``position``
    (relative to the star on ``start_date``), one sample per day."""
    if ndays is None:
        ndays = int(round(365.25 * (nyears or 1.0)))
    t0 = Time(start_date)
    times = t0 + np.arange(ndays) * u.day
    app = _apparent(star, times)
    ref = _apparent(star, t0)
    # background object stays at the reference apparent position; offset = bg − star
    d_ra = (ref.ra - app.ra).to(u.arcsec).value * np.cos(ref.dec.rad)
    d_dec = (ref.dec - app.dec).to(u.arcsec).value
    d_e = position[0] + d_ra
    d_n = position[1] + d_dec
    tr = BackgroundTrack(times, d_e, d_n)
    if posfile is not None:
        e, n, se, sn, dates = read_position_file(posfile)
        te = Time(dates)
        appe = _apparent(star, te)
        ref0 = _apparent(star, te[0])
        pe = e[0] + (ref0.ra - appe.ra).to(u.arcsec).value * np.cos(ref0.dec.rad)
        pn = n[0] + (ref0.dec - appe.dec).to(u.arcsec).value
        tr.predicted_at_epochs = np.column_stack([pe, pn])
        if e.size > 1:
            tr.chisq = float(
                np.sum(((e[1:] - pe[1:]) / se[1:]) ** 2 + ((n[1:] - pn[1:]) / sn[1:]) ** 2)
            )
            log.info(
                "background hypothesis chi^2 = %.2f over %d later epochs", tr.chisq, e.size - 1
            )
    if plot is not None:
        from charis_dpp.plotting import finish, plt

        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.plot(d_e, d_n, "-", color="0.3", label="background object track")
        if posfile is not None:
            ax.errorbar(e[0], n[0], sn[0], se[0], fmt="o", color="C3", label="first epoch")
            if e.size > 1:
                ax.errorbar(e[1:], n[1:], sn[1:], se[1:], fmt="o", color="C2", label="later epochs")
                ax.plot(pe[1:], pn[1:], "x", color="k", label="predicted if background")
        ax.invert_xaxis()
        ax.set_xlabel("E offset [arcsec]")
        ax.set_ylabel("N offset [arcsec]")
        ax.set_aspect("equal", adjustable="datalim")
        ax.legend(frameon=False, fontsize=8)
        finish(fig, Path(plot))
    return tr
