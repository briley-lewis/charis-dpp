"""WCS bookkeeping: north-up rotation and the post-NBS parity correction.

``northup`` replaces ``charis_northup.pro`` + ``hrot2.pro``: the data are
rotated once by ``rotation_from_header(h) - north_offset`` and the CD matrix /
CRPIX are updated once (the IDL version rotated the header twice, BUGS.md #269).
"""

from __future__ import annotations

import numpy as np
from astropy.io.fits import Header

from charis_dpp.geometry.grids import image_center
from charis_dpp.image.interp import rotate_cube, rotate_image

__all__ = ["cd_matrix", "rotation_from_header", "rotate_wcs", "northup", "nbs_parity_correct"]


def cd_matrix(h: Header) -> np.ndarray | None:
    """CD matrix ``[[CD1_1, CD1_2], [CD2_1, CD2_2]]`` from CD or CDELT/CROTA cards."""
    if all(k in h for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2")):
        return np.array([[h["CD1_1"], h["CD1_2"]], [h["CD2_1"], h["CD2_2"]]], dtype=float)
    if "CDELT1" in h and "CDELT2" in h:
        rot = np.deg2rad(float(h.get("CROTA2", 0.0)))
        c1, c2 = float(h["CDELT1"]), float(h["CDELT2"])
        return np.array(
            [[c1 * np.cos(rot), -c2 * np.sin(rot)], [c1 * np.sin(rot), c2 * np.cos(rot)]]
        )
    return None


def rotation_from_header(h: Header) -> float:
    """Astrolib ``getrot``: rotation of the image from north-up, degrees (port of ``getrot.pro``).

    ``extast`` stores ``astr.cd[i, j] = CD(i+1)_(j+1)``; ``sgn`` is the sign of
    the CD determinant (−1 for the parity-flipped CHARIS frames),
    ``rot = atan(sgn·CD1_2, sgn·CD1_1)`` (x axis), ``rot2 = atan(−CD2_1, CD2_2)``
    (y axis), averaged when they differ by less than 2 deg, else ``rot`` is kept.
    Checked against IDL ``getrot`` on unflipped, flipped and asymmetric headers.
    The earlier ``atan2(CD1_2, CD2_2)`` ignored ``sgn`` and so had the opposite
    sign for CHARIS (``+3.31`` for ``−3.31`` deg): every north-up product was
    rotated by twice the field angle (BUGS.md #332).
    """
    cd = cd_matrix(h)
    if cd is None:
        if "TOT_ROT" in h:
            return float(h["TOT_ROT"])
        raise KeyError("header has no CD matrix / CDELT to derive the rotation from")
    cd1_1, cd1_2 = float(cd[0, 0]), float(cd[0, 1])
    cd2_1, cd2_2 = float(cd[1, 0]), float(cd[1, 1])
    det = cd1_1 * cd2_2 - cd2_1 * cd1_2
    sgn = -1.0 if det < 0 else 1.0
    if cd1_2 == 0 and cd2_1 == 0:
        return 0.0
    rot = np.arctan2(sgn * cd1_2, sgn * cd1_1)
    rot2 = np.arctan2(-cd2_1, cd2_2)
    if rot != rot2:
        tpi = 2 * np.pi
        if abs(rot - rot2) < np.deg2rad(2):
            rot = 0.5 * (rot + rot2)
        elif abs(rot - rot2 - tpi) < np.deg2rad(2):
            rot = 0.5 * (rot + rot2 - tpi)
        elif abs(rot - rot2 + tpi) < np.deg2rad(2):
            rot = 0.5 * (rot + rot2 + tpi)
    return float(np.degrees(rot))


def rotate_wcs(h: Header, deg: float, center: tuple[float, float]) -> None:
    """Update CD/CRPIX in place for a data rotation by ``deg`` about ``center`` (0-based), ``rotat`` conventions."""
    cd = cd_matrix(h)
    if cd is None:
        h["ROTINFO"] = (deg, "data rotated; no WCS to update")
        return
    t = np.deg2rad(deg)
    ct, st = np.cos(t), np.sin(t)
    # rotat forward map: pf - c = R (pi - c), R = [[ct, st], [-st, ct]]
    R = np.array([[ct, st], [-st, ct]])
    Rinv = R.T
    cd_new = cd @ Rinv
    x0, y0 = center
    crpix = np.array([float(h.get("CRPIX1", x0 + 1)), float(h.get("CRPIX2", y0 + 1))]) - 1.0
    crpix_new = np.array([x0, y0]) + R @ (crpix - np.array([x0, y0]))
    h["CD1_1"], h["CD1_2"] = cd_new[0, 0], cd_new[0, 1]
    h["CD2_1"], h["CD2_2"] = cd_new[1, 0], cd_new[1, 1]
    h["CRPIX1"], h["CRPIX2"] = crpix_new[0] + 1.0, crpix_new[1] + 1.0
    for k in ("CDELT1", "CDELT2", "CROTA1", "CROTA2"):
        h.pop(k, None)


def northup(
    data: np.ndarray,
    h1: Header,
    north_offset_deg: float,
    *,
    center: tuple[float, float] | None = None,
    frame_index: int | None = None,
) -> tuple[np.ndarray, float]:
    """Rotate a frame or cube north-up.  Returns ``(rotated, rotang)`` and updates ``h1``.

    ``rotang`` is the angle recorded as ``ROTANG`` (and ``ROTANG_n`` when
    ``frame_index`` is given): the rotation applied, in ``rotat`` convention
    (``-(npa) + north_offset`` → data rotated by ``npa - north_offset``).
    """
    npa = rotation_from_header(h1)
    d_par_ang = -npa + north_offset_deg
    rot_deg = -d_par_ang
    xc, yc = image_center(data.shape) if center is None else center
    out = (
        rotate_cube(data, rot_deg, (xc, yc))
        if data.ndim == 3
        else rotate_image(data, rot_deg, (xc, yc))
    )
    rotate_wcs(h1, rot_deg, (xc, yc))
    h1["ROTANG"] = (rot_deg, "rotation applied for north-up [deg]")
    h1["NORTHOFF"] = (north_offset_deg, "north offset included in ROTANG [deg]")
    if frame_index is not None:
        h1[f"ROTANG_{frame_index}"] = (rot_deg, "per-frame north-up rotation [deg]")
    return out, rot_deg


def nbs_parity_correct(
    cube: np.ndarray, h0: Header, h1: Header, rotzp_deg: float, parang: float | None = None
) -> np.ndarray:
    """Post-NBS parity flip + WCS rebuild (``charis_nbs_corr.pro``): x-flip, CD from ``-(PARANG + 113 - rotzp)``."""
    extrarot = 113.0 - rotzp_deg
    if parang is None:
        parang = float(h1.get("PARANG", h0.get("PARANG", 0.0)))
    totrot = -1.0 * (parang + extrarot)
    xp = float(h1.get("XPIXSCAL", h0.get("XPIXSCAL", 0.0)))
    yp = float(h1.get("YPIXSCAL", h0.get("YPIXSCAL", 0.0)))
    a = np.deg2rad(totrot)
    for h in (h1, h0):
        for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2", "XPIXSCAL", "YPIXSCAL", "TOT_ROT"):
            h.pop(k, None)
        h["PARITY"] = (1, "flipped from N-u,E-r to N-u, E-l (post-NBS)")
        h["CD1_1"] = np.cos(a) * xp
        h["CD1_2"] = np.sin(a) * xp
        h["CD2_1"] = -np.sin(a) * yp
        h["CD2_2"] = np.cos(a) * yp
        h["XPIXSCAL"] = xp
        h["YPIXSCAL"] = yp
        h["ROT_ZP"] = (rotzp_deg, "zp rot (degrees)")
        h["EX_ROT"] = (extrarot, "extrarot:113-zp (degrees)")
        h["TOT_ROT"] = (totrot, "Total rotation angle (degrees)")
        h.add_comment("post-NBS Parity Flip")
    return cube[..., ::-1]
