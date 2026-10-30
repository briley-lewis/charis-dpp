"""Mueller-matrix model, Stokes utilities and HWP matching."""

from __future__ import annotations

import numpy as np
from charis_dpp.config.instrument import instrument_for
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.pdi.hwp_match import HwpTable, assign_cycles, read_hwp_log
from charis_dpp.pdi.mueller.components import (
    derotator_retardances,
    hwp_retardance,
    m3_model,
    wollaston_efficiency,
)
from charis_dpp.pdi.mueller.matrix import (
    charis_mueller_matrix,
    m_component,
    m_elliptical_retarder,
    m_rotation,
    m_wollaston,
)
from charis_dpp.pdi.stokes import (
    adjust_pol_axis,
    average_angle,
    azimuthal_stokes,
    polarisation_angle_grid,
)


def test_rotation_and_ideal_components():
    assert np.allclose(m_rotation(37.0) @ m_rotation(-37.0), np.eye(4))
    assert np.allclose(m_component(0.0, 0.0), np.eye(4))
    assert np.allclose(m_elliptical_retarder(0.0, 0.0, 0.0, 0.0), np.eye(4))
    # ideal half-wave plate flips U and V
    assert np.allclose(m_component(0.0, 180.0), np.diag([1, 1, -1, -1]))
    # ideal Wollaston: right - left = full Q analyser row
    diff = m_wollaston(1.0, 1.0) - m_wollaston(1.0, -1.0)
    assert np.allclose(diff[0], [0, 1, 0, 0])
    # ideal train with HWP at 0 and 45 deg: Q and -Q analysers
    row0 = charis_mueller_matrix(1.0, 0, 0, 0, 0, 0, 180.0, 0, 0, 0, 0.0, 0.0, 0.0)[0]
    row45 = charis_mueller_matrix(1.0, 0, 0, 0, 0, 0, 180.0, 0, 0, 0, 45.0, 0.0, 0.0)[0]
    row22 = charis_mueller_matrix(1.0, 0, 0, 0, 0, 0, 180.0, 0, 0, 0, 22.5, 0.0, 0.0)[0]
    assert np.allclose(row0, [0, 1, 0, 0]) and np.allclose(row45, [0, -1, 0, 0])
    assert np.allclose(np.abs(row22), [0, 0, 1, 0])


def test_component_models_epochs():
    lam = np.linspace(1.16, 2.37, 22)
    pre = instrument_for(58000.0)
    post = instrument_for("2026-06-01")
    d_hwp, eps = hwp_retardance(lam, pre)
    assert eps == 0.0 and np.all(np.abs(d_hwp - 180.0) < 25.0)  # near-achromatic HWP
    d_pre, e_pre = m3_model(lam, pre)
    d_post, e_post = m3_model(lam, post)
    assert np.all((d_pre > 150) & (d_pre < 200)) and np.all(np.abs(e_pre) < 0.1)
    assert not np.allclose(d_pre, d_post)  # the recoat epoch switches the fit
    dh, d45, dr, e = derotator_retardances(lam, pre)
    assert dh.shape == lam.shape and np.all(dh > 0) and e == 0.0
    w = wollaston_efficiency(lam, pre)
    assert np.all((w > 0.9) & (w <= 1.0))


def test_azimuthal_stokes_roundtrip():
    shape = (61, 61)
    phi = polarisation_angle_grid(shape)
    r = radial_grid(shape)
    ring = np.exp(-0.5 * ((r - 15) / 2.0) ** 2)
    # tangential polarisation: Qphi = ring, Uphi = 0 -> Q = -Qphi cos2phi, U = -Qphi sin2phi
    q = -ring * np.cos(2 * phi)
    u = -ring * np.sin(2 * phi)
    qphi, uphi = azimuthal_stokes(q, u, phi)
    assert np.allclose(qphi, ring, atol=1e-12) and np.allclose(uphi, 0.0, atol=1e-12)
    # rotating the polarisation axis by 90 deg flips the sign of Q and U
    q2, u2 = adjust_pol_axis(q, u, 90.0, phi)
    assert np.allclose(q2, -q, atol=1e-12) and np.allclose(u2, -u, atol=1e-12)
    q0, u0 = adjust_pol_axis(q, u, 0.0, phi)
    assert np.allclose(q0, q) and np.allclose(u0, u)
    assert (
        abs(average_angle([350.0, 10.0])) < 1e-9 or abs(average_angle([350.0, 10.0]) - 360) < 1e-9
    )


def test_hwp_log_and_cycles(tmp_path):
    log = tmp_path / "vampHWPLog_20201003.csv"
    log.write_text(
        "angle,ut_start,ut_end,target\n"
        "0.0,2020-10-03T10:00:00,2020-10-03T10:00:03,AB Aur\n"
        "22.5,10:01:00,10:01:03,AB Aur\n"
        "45,10:02:00,,AB Aur\n"
        "67.5,2020-10-03T100300,2020-10-03T100303,AB Aur\n"
    )
    angs, t0, t1 = read_hwp_log(log)
    assert np.allclose(angs, [0, 22.5, 45, 67.5])
    assert np.allclose(t0, [10.0, 10 + 1 / 60, 10 + 2 / 60, 10 + 3 / 60])
    assert np.isnan(t1[2])
    # 9 frames: two complete cycles plus a stray 0-deg frame
    pos = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0])
    times = np.arange(9) / 60.0
    cyc = assign_cycles(pos, times, 4)
    assert list(cyc) == [0, 0, 0, 0, 1, 1, 1, 1, -1]
    t = HwpTable([f"n{i:04d}leftreg.fits" for i in range(9)], list(range(9)), pos, pos * 22.5, cyc)
    p = t.save(tmp_path / "hwp_info.txt")
    t2 = HwpTable.load(p)
    assert t2.ncycles == 2 and list(t2.hwp_cycle) == list(cyc) and t2.filenums == list(range(9))
