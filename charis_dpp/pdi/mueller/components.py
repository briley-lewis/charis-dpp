"""Component models (``pdi/mueller_matrix/component_models/*.pro``).

Every parameter comes from the instrument epoch table (``[epoch.pdi.mueller.*]``),
so the pre- and post-recoat M3 fits are selected by date (BUGS.md #307) and the
HWP birefringence is always evaluated from the Sellmeier equations instead of a
22-entry precomputed table (#309).  Wavelengths are in µm throughout; angles
returned in degrees unless stated.
"""

from __future__ import annotations

import numpy as np

from charis_dpp.config.instrument import Instrument
from charis_dpp.config.paths import package_data
from charis_dpp.io.tables import read_csv_columns

__all__ = [
    "sellmeier_ghosh",
    "sellmeier_terms",
    "hwp_retardance",
    "complex_refractive_index",
    "fresnel_rs_rp",
    "m3_model",
    "derotator_retardances",
    "wollaston_efficiency",
]


def sellmeier_ghosh(lam_um, coeff) -> np.ndarray:
    """``n² = A + B·L²/(L²−C) + D·L²/(L²−E)`` (Ghosh 1999 form; ``refractive_index.pro``)."""
    a, b, c, d, e = (float(v) for v in coeff)
    l2 = np.asarray(lam_um, dtype=float) ** 2
    return np.sqrt(a + b * l2 / (l2 - c) + d * l2 / (l2 - e))


def sellmeier_terms(lam_um, terms) -> np.ndarray:
    """``n² = 1 + Σ B_i·L²/(L² − C_i²)`` (refractiveindex.info form, Dodge 1984 for MgF2)."""
    l2 = np.asarray(lam_um, dtype=float) ** 2
    n2 = np.ones_like(l2)
    for b, c in terms:
        n2 = n2 + float(b) * l2 / (l2 - float(c) ** 2)
    return np.sqrt(n2)


def hwp_retardance(lam_um, inst: Instrument) -> tuple[np.ndarray, float]:
    """(retardance δ_HWP in degrees per wavelength, diattenuation ε_HWP) of the SiO2+MgF2 achromatic HWP.

    Hart et al. 2021 eq. 23: δ = 2π/λ · (W_MgF2·Δn_MgF2 − W_SiO2·Δn_SiO2).
    """
    h = inst.get("pdi.mueller.hwp")
    lam = np.asarray(lam_um, dtype=float)
    no_sio2 = sellmeier_ghosh(lam, h.sio2_o)
    ne_sio2 = sellmeier_ghosh(lam, h.sio2_e)
    no_mgf2 = sellmeier_terms(lam, h.mgf2_o)
    ne_mgf2 = sellmeier_terms(lam, h.mgf2_e)
    d_sio2 = float(h.w_sio2_mm) * (ne_sio2 - no_sio2)
    d_mgf2 = float(h.w_mgf2_mm) * (ne_mgf2 - no_mgf2)
    lam_mm = lam * 1e-3
    delta = (2 * np.pi / lam_mm) * (d_mgf2 - d_sio2)
    return np.degrees(delta), float(h.get("epsilon", 0.0))


def complex_refractive_index(eps_r1, eps_r2) -> np.ndarray:
    """``n + ik`` from the real and imaginary relative permittivity."""
    eps_r1 = np.asarray(eps_r1, dtype=float)
    eps_r2 = np.asarray(eps_r2, dtype=float)
    mag = np.sqrt(eps_r1**2 + eps_r2**2)
    n = np.sqrt(mag + eps_r1) / np.sqrt(2.0)
    k = np.sqrt(mag - eps_r1) / np.sqrt(2.0)
    return n + 1j * k


def fresnel_rs_rp(theta_i_rad: float, n1, n2) -> tuple[np.ndarray, np.ndarray]:
    """Fresnel amplitude reflection coefficients (s, p) at incidence ``theta_i`` from ``n1`` into ``n2``."""
    cos_i = np.cos(theta_i_rad)
    cos_t = np.sqrt(1.0 - ((n1 / n2) * np.sin(theta_i_rad)) ** 2 + 0j)
    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return r_s, r_p


def m3_model(lam_um, inst: Instrument) -> tuple[np.ndarray, np.ndarray]:
    """(retardance δ_M3 in degrees, diattenuation ε_M3) of the tertiary mirror from the
    epoch's complex-permittivity power-law fit (``m3_model.pro``)."""
    m = inst.get("pdi.mueller.m3")
    lam_m = np.asarray(lam_um, dtype=float) * 1e-6
    eps_r1 = -(10.0 ** float(m.b1)) * lam_m ** float(m.m1)
    eps_r2 = (10.0 ** float(m.b2)) * lam_m ** float(m.m2)
    n2 = complex_refractive_index(eps_r1, eps_r2)
    n1 = float(m.get("n_air", 1.0))
    r_s, r_p = fresnel_rs_rp(np.deg2rad(float(m.get("incidence_deg", 45.0))), n1, n2)
    rs2, rp2 = np.abs(r_s) ** 2, np.abs(r_p) ** 2
    eps = (rs2 - rp2) / (rs2 + rp2)
    delta = np.mod(np.angle(r_s) - np.angle(r_p) + 2 * np.pi, 2 * np.pi)
    return np.degrees(delta), eps


def _wavelength_fits(inst: Instrument) -> dict[str, np.ndarray]:
    rel = str(inst.get("pdi.mueller.wavelength_fits_csv"))
    cols = read_csv_columns(package_data(rel))
    return {k.lower(): np.asarray(v, dtype=float) for k, v in cols.items() if k.strip()}


def derotator_retardances(
    lam_um, inst: Instrument
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Elliptical-retarder parameters of the image rotator ``(δ_h, δ_45, δ_r)`` in radians
    (as fitted by McIntosh+2026, interpolated in wavelength) and its diattenuation."""
    t = _wavelength_fits(inst)
    lam_nm = np.asarray(lam_um, dtype=float) * 1e3
    w = t["wavelength_bin"]
    dh = np.interp(lam_nm, w, t["image_rotator_phi_h"])
    d45 = np.interp(lam_nm, w, t["image_rotator_phi_45"])
    dr = np.interp(lam_nm, w, t["image_rotator_phi_r"])
    return dh, d45, dr, float(inst.get("pdi.mueller.derotator.epsilon", 0.0))


def wollaston_efficiency(lam_um, inst: Instrument) -> np.ndarray:
    """Wollaston diattenuation ε_bs per wavelength (McIntosh+2026 fit, interpolated)."""
    t = _wavelength_fits(inst)
    lam_nm = np.asarray(lam_um, dtype=float) * 1e3
    return np.interp(lam_nm, t["wavelength_bin"], t["wollaston_eta"])
