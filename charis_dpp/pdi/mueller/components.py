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


def interpol_idl(x, xp, fp) -> np.ndarray:
    """IDL ``interpol(fp, xp, x)``: piecewise linear *with* linear extrapolation
    beyond both ends (``np.interp`` clamps to the end values; the CHARIS
    broadband channels 0 and 21 lie 0.4 nm outside the McIntosh+2026 table)."""
    x = np.asarray(x, dtype=float)
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    y = np.interp(x, xp, fp)
    lo = x < xp[0]
    if lo.any():
        y[lo] = fp[0] + (x[lo] - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0])
    hi = x > xp[-1]
    if hi.any():
        y[hi] = fp[-1] + (x[hi] - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
    return y


# hwp_model.pro uses these precomputed indices for the 22-channel broadband
# case instead of evaluating the Sellmeier equations (they differ from the
# formulas in the 4th decimal, 0.04 deg of HWP retardance); mirrored here so
# the two sides agree (BUGS.md #331).
_HWP_BROADBAND_NO_MGF2 = np.array([1.37264460186324, 1.37242217890004, 1.3721950699442604, 1.3719621982496155, 1.3717224594867308, 1.3714747161454446, 1.3712177917570674, 1.3709504649005464, 1.3706714629536383, 1.370379455549121, 1.3700730476930256, 1.3697507724980869, 1.3694110834841233, 1.3690523463906459, 1.3686728304441675, 1.368270699015774, 1.3678439995997282, 1.3673906530353213, 1.366908441887303, 1.3663949978894672, 1.3658477883469353, 1.3652641013782865])  # fmt: skip
_HWP_BROADBAND_NE_MGF2 = np.array([1.384198687268396, 1.3839628415351612, 1.3837218336935047, 1.383474521071655, 1.3832197307469136, 1.382956253598398, 1.3826828381620595, 1.382398184249572, 1.382100936289585, 1.3817896763486257, 1.3814629167856005, 1.381119092489746, 1.3807565526501837, 1.380373551998335, 1.3799682414613015, 1.37953865815685, 1.3790827146553932, 1.3785981874250786, 1.3780827043685924, 1.3775337313486054, 1.3769485575889246, 1.3763242798227924])  # fmt: skip
_HWP_BROADBAND_NO_SIO2 = np.array([1.532817938541442, 1.5322934166189457, 1.531754934215596, 1.531199788598694, 1.5306251760450662, 1.5300281741782964, 1.5294057231440228, 1.5287546054184482, 1.5280714240176634, 1.527352578848145, 1.5265942409020412, 1.525792323956521, 1.524942453390836, 1.524039931668601, 1.5230796999650995, 1.5220562953290564, 1.5209638026676928, 1.5197958007146606, 1.5185453009917602, 1.5172046785872932, 1.5157655933534087, 1.51421889984422])  # fmt: skip
_HWP_BROADBAND_NE_SIO2 = np.array([1.5415034516622548, 1.5409579239812228, 1.5403976973079714, 1.539819958217567, 1.5392217872099536, 1.5386001402848022, 1.537951829297335, 1.5372735008806506, 1.536561613691175, 1.5358124137052376, 1.5350219072561897, 1.5341858314549301, 1.5332996215887311, 1.5323583750238043, 1.531356811065872, 1.5302892261380394, 1.529149443529453, 1.527930756832346, 1.5266258660286436, 1.5252268049894229, 1.5237248589184071, 1.5221104699752577])  # fmt: skip


def hwp_retardance(lam_um, inst: Instrument) -> tuple[np.ndarray, float]:
    """(retardance δ_HWP in degrees per wavelength, diattenuation ε_HWP) of the SiO2+MgF2 achromatic HWP.

    Hart et al. 2021 eq. 23: δ = 2π/λ · (W_MgF2·Δn_MgF2 − W_SiO2·Δn_SiO2).
    """
    h = inst.get("pdi.mueller.hwp")
    lam = np.asarray(lam_um, dtype=float)
    if lam.size == 22:  # hwp_model.pro: precomputed broadband indices
        no_sio2, ne_sio2 = _HWP_BROADBAND_NO_SIO2, _HWP_BROADBAND_NE_SIO2
        no_mgf2, ne_mgf2 = _HWP_BROADBAND_NO_MGF2, _HWP_BROADBAND_NE_MGF2
    else:
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
    dh = interpol_idl(lam_nm, w, t["image_rotator_phi_h"])
    d45 = interpol_idl(lam_nm, w, t["image_rotator_phi_45"])
    dr = interpol_idl(lam_nm, w, t["image_rotator_phi_r"])
    return dh, d45, dr, float(inst.get("pdi.mueller.derotator.epsilon", 0.0))


def wollaston_efficiency(lam_um, inst: Instrument) -> np.ndarray:
    """Wollaston diattenuation ε_bs per wavelength (McIntosh+2026 fit, interpolated)."""
    t = _wavelength_fits(inst)
    lam_nm = np.asarray(lam_um, dtype=float) * 1e3
    return interpol_idl(lam_nm, t["wavelength_bin"], t["wollaston_eta"])
