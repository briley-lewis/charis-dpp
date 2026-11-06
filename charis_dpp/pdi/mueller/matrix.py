"""Mueller matrices (``pdi/mueller_matrix/matrix_terms/*.pro``, ``charis_pdi_mueller_matrix.pro``,
``get_charis_mueller_matrix.pro``).

Matrices are numpy ``(4, 4)`` arrays indexed ``[row, col]`` acting on
``(I, Q, U, V)``; the IDL ``##`` products map to ``@`` in the same order.
Angles in degrees except the elliptical-retarder retardances (radians, as
fitted).

Fixes vs. IDL: the derotator angle and parallactic angle are read from the
pipeline header (``imprep`` copies ``D_IMRANG`` from the raw cube) instead of
re-opening ``data/raw`` (#318); the site latitude comes from the instrument
table (#306).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.io.fits import Header

from charis_dpp.config.instrument import Instrument
from charis_dpp.io.fits import header_get
from charis_dpp.pdi.mueller.components import (
    derotator_retardances,
    hwp_retardance,
    m3_model,
    wollaston_efficiency,
)

__all__ = [
    "m_wollaston",
    "m_component",
    "m_elliptical_retarder",
    "m_rotation",
    "charis_mueller_matrix",
    "hwp_sync_angle",
    "mueller_rows_per_channel",
    "MuellerTerms",
]


def m_wollaston(eps: float, half: float) -> np.ndarray:
    """One Wollaston beam (``half = ±1``) with diattenuation ``eps``."""
    m = np.zeros((4, 4))
    m[0, 0] = m[1, 1] = 0.5
    m[0, 1] = m[1, 0] = 0.5 * half * eps
    m[2, 2] = m[3, 3] = 0.5 * np.sqrt(1.0 - eps**2)
    return m


def m_component(eps: float, delta_deg: float) -> np.ndarray:
    """Linear retarder/diattenuator with diattenuation ``eps`` and retardance ``delta`` (deg)."""
    d = np.deg2rad(delta_deg)
    s = np.sqrt(1.0 - eps**2)
    m = np.zeros((4, 4))
    m[0, 0] = m[1, 1] = 1.0
    m[0, 1] = m[1, 0] = eps
    m[2, 2] = m[3, 3] = s * np.cos(d)
    m[2, 3] = s * np.sin(d)
    m[3, 2] = -s * np.sin(d)
    return m


def m_elliptical_retarder(eps: float, dh: float, d45: float, dr: float) -> np.ndarray:
    """General elliptical retarder with retardance vector ``(δ_h, δ_45, δ_r)`` in radians."""
    phi2 = dh**2 + dr**2 + d45**2
    m = np.zeros((4, 4))
    m[0, 0] = 1.0
    if phi2 == 0:
        m[1, 1] = m[2, 2] = m[3, 3] = 1.0
        return m
    phi = np.sqrt(phi2)
    c, s = np.cos(phi), np.sin(phi)
    t = 1.0 - c
    m[1, 1] = (dh**2 + (d45**2 + dr**2) * c) / phi2
    m[1, 2] = (d45 * dh * t) / phi2 + (dr * s) / phi
    m[1, 3] = (dr * dh * t) / phi2 - (d45 * s) / phi
    m[2, 1] = (dh * d45 * t) / phi2 - (dr * s) / phi
    m[2, 2] = (d45**2 + (dh**2 + dr**2) * c) / phi2
    m[2, 3] = (dr * d45 * t) / phi2 + (dh * s) / phi
    m[3, 1] = (dh * dr * t) / phi2 + (d45 * s) / phi
    m[3, 2] = (d45 * dr * t) / phi2 - (dh * s) / phi
    m[3, 3] = (dr**2 + (dh**2 + d45**2) * c) / phi2
    return m


def m_rotation(theta_deg: float) -> np.ndarray:
    """Rotation of the (Q, U) frame by ``theta`` (deg)."""
    t = 2.0 * np.deg2rad(theta_deg)
    m = np.eye(4)
    m[1, 1] = m[2, 2] = np.cos(t)
    m[1, 2] = np.sin(t)
    m[2, 1] = -np.sin(t)
    return m


def charis_mueller_matrix(
    eps_bs: float,
    eps_der: float,
    dh_der: float,
    d45_der: float,
    dr_der: float,
    eps_hwp: float,
    delta_hwp: float,
    eps_tel: float,
    delta_tel: float,
    theta_der: float,
    theta_hwp: float,
    altitude: float,
    parang: float,
    *,
    parity_flip: bool = False,
) -> np.ndarray:
    """Double-difference Mueller matrix ``M_right − M_left`` of telescope + HWP + derotator + Wollaston.

    ``parity_flip`` inserts ``diag(1, 1, −1, −1)`` after the telescope term (an
    x-flip of the image reverses the sign of U and V); the fork omits it
    because the NBS parity correction has already flipped post-2025 data
    (BUGS.md #308).
    """
    r_der = (
        m_rotation(-theta_der)
        @ m_elliptical_retarder(eps_der, dh_der, d45_der, dr_der)
        @ m_rotation(theta_der)
    )
    r_hwp = m_rotation(-theta_hwp) @ m_component(eps_hwp, delta_hwp) @ m_rotation(theta_hwp)
    r_tel = m_rotation(-altitude) @ m_component(eps_tel, delta_tel) @ m_rotation(parang)
    if parity_flip:
        r_tel = np.diag([1.0, 1.0, -1.0, -1.0]) @ r_tel
    m0 = r_der @ r_hwp @ r_tel
    return m_wollaston(eps_bs, 1.0) @ m0 - m_wollaston(eps_bs, -1.0) @ m0


def hwp_sync_angle(altitude: float, azimuth: float, lat_deg: float) -> float:
    """HWP angle offset ``0.5·field_rotation + altitude`` used in ADI-synchronised mode."""
    a_tel = np.deg2rad(lat_deg)
    az = np.deg2rad(azimuth - 180.0)
    alt = np.deg2rad(altitude)
    fr = np.degrees(np.arctan2(np.sin(az), np.sin(alt) * np.cos(az) + np.cos(alt) * np.tan(a_tel)))
    return 0.5 * fr + altitude


class MuellerTerms:
    """Per-wavelength component parameters of one epoch, computed once."""

    def __init__(self, wav_um: np.ndarray, inst: Instrument):
        self.wav = np.asarray(wav_um, dtype=float)
        self.inst = inst
        self.delta_hwp, self.eps_hwp = hwp_retardance(self.wav, inst)
        self.dh_der, self.d45_der, self.dr_der, self.eps_der = derotator_retardances(self.wav, inst)
        self.delta_tel, self.eps_tel = m3_model(self.wav, inst)
        self.eps_bs = wollaston_efficiency(self.wav, inst)
        mm = inst.get("pdi.mueller")
        self.dtheta_hwp = float(mm.get("dtheta_hwp_deg", 0.0))
        self.dtheta_der = float(mm.get("dtheta_der_deg", 0.0))
        self.parity_flip = bool(mm.get("flip_matrix", False))


def mueller_rows_per_channel(
    h: Header,
    hwpang: float,
    terms: MuellerTerms,
    *,
    adisync: bool = True,
    h0: Header | None = None,
    raw_dir: Path | None = None,
) -> np.ndarray:
    """First rows ``(nλ, 4)`` of the double-difference Mueller matrix for one frame.

    ``h`` is the pipeline extension header carrying ``ALTITUDE``, ``AZIMUTH``
    and ``PARANG``; ``D_IMRANG`` (the image-rotator angle, which the IDL reads
    from the raw cube's third extension) is looked up in ``h`` and then in the
    primary header ``h0``, where ``imprep`` copies it, and finally in the raw
    cube's third extension under ``raw_dir`` (what the IDL does, for frames
    prepared by the IDL ``imprep``).  A missing ``D_IMRANG`` is an error: silently using 0 deg gave a wrong correction (BUGS.md #330).
    ``hwpang`` is the HWP angle from the HWP table.
    """
    altitude = float(header_get(h, "ALTITUDE", 0.0))
    parang = float(header_get(h, "PARANG", header_get(h, "PA", 0.0)))
    theta_der = header_get(h, "D_IMRANG", None)
    if theta_der is None and h0 is not None:
        theta_der = header_get(h0, "D_IMRANG", None)
    if theta_der is None and raw_dir is not None:
        # frames prepared by the IDL imprep carry no D_IMRANG: do what
        # get_charis_mueller_matrix.pro does and read the raw cube's 3rd extension
        orig = header_get(h, "ORIGNAME", None)
        raw = Path(raw_dir) / f"CRSA{orig}_cube.fits" if orig is not None else None
        if raw is not None and raw.exists():
            theta_der = header_get(fits.getheader(raw, 3), "D_IMRANG", None)
    if theta_der is None:
        raise KeyError("D_IMRANG (image rotator angle) missing from the frame headers")
    theta_der = float(theta_der)
    if adisync:
        azimuth = float(header_get(h, "AZIMUTH", 180.0))
        theta_hwp = hwp_sync_angle(altitude, azimuth, terms.inst.latitude_deg) + hwpang
    else:
        theta_hwp = hwpang
    theta_der += terms.dtheta_der
    theta_hwp += terms.dtheta_hwp
    rows = np.zeros((terms.wav.size, 4))
    for il in range(terms.wav.size):
        m = charis_mueller_matrix(
            float(terms.eps_bs[il]),
            terms.eps_der,
            float(terms.dh_der[il]),
            float(terms.d45_der[il]),
            float(terms.dr_der[il]),
            terms.eps_hwp,
            float(terms.delta_hwp[il]),
            float(terms.eps_tel[il]),
            float(terms.delta_tel[il]),
            theta_der,
            theta_hwp,
            altitude,
            parang,
            parity_flip=terms.parity_flip,
        )
        rows[il] = m[0, :]
    return rows
