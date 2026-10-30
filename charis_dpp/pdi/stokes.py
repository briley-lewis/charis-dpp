"""Stokes utilities (``pdi/utilities/*.pro``)."""

from __future__ import annotations

import numpy as np

__all__ = ["azimuthal_stokes", "adjust_pol_axis", "average_angle", "polarisation_angle_grid"]


def polarisation_angle_grid(shape, center=None) -> np.ndarray:
    """Azimuth φ used by ``charis_pdi`` (``reverse(transpose(angarr(n)))``): the angle
    of each pixel measured from the +y axis (north) towards −x (east on a
    north-up, east-left image), in radians."""
    from charis_dpp.geometry.grids import xy_grids

    dx, dy = xy_grids(shape, center)
    return np.arctan2(-dx, dy)


def azimuthal_stokes(
    q: np.ndarray, u: np.ndarray, phi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """``(Qφ, Uφ)`` with ``Qφ = −Q cos2φ − U sin2φ``, ``Uφ = Q sin2φ − U cos2φ`` (Schmid+2006 convention as in the IDL)."""
    c, s = np.cos(2.0 * phi), np.sin(2.0 * phi)
    qphi = -q * c - u * s
    uphi = q * s - u * c
    return qphi, uphi


def adjust_pol_axis(
    q: np.ndarray, u: np.ndarray, da_deg: float, phi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate the polarisation reference axis by ``da`` degrees (``adjust_pol_axis.pro``):
    the azimuthal Stokes are re-projected with ``φ + da``."""
    qphi, uphi = azimuthal_stokes(q, u, phi)
    p = phi + np.deg2rad(da_deg)
    c, s = np.cos(2.0 * p), np.sin(2.0 * p)
    q_out = uphi * s - qphi * c
    u_out = -uphi * c - qphi * s
    return q_out, u_out


def average_angle(angles_deg, axis=None) -> float | np.ndarray:
    """Circular mean of angles in degrees, in ``[0, 360)``."""
    a = np.deg2rad(np.asarray(angles_deg, dtype=float))
    x = np.mean(np.cos(a), axis=axis)
    y = np.mean(np.sin(a), axis=axis)
    return np.mod(np.degrees(np.arctan2(y, x)) + 360.0, 360.0)
