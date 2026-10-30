"""Unit tests for the forward-modelling helpers, chi^2 utilities and the covariance model."""

from __future__ import annotations

import numpy as np
import pytest
from charis_dpp.partialsub.diskfit import chisq_map, model_in_frame, roi_mask, scale_model_ls
from charis_dpp.partialsub.inject import frame_position
from charis_dpp.specphotcal.spec_covar import covariance_matrix, gb_model
from charis_dpp.tools.chisq import binned_chisq, chisq, chisq_covar


def test_model_in_frame_matches_frame_position():
    shape = (81, 81)
    center = (40.0, 40.0)
    model = np.zeros(shape)
    x0, y0 = 55.0, 48.0
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    model += np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / 1.5**2)
    for rot in (0.0, 37.0, -112.5, 90.0):
        im = model_in_frame(model, rot, center)
        iy, ix = np.unravel_index(np.nanargmax(im), im.shape)
        px, py = frame_position(x0, y0, rot, center)
        assert abs(ix - px) <= 1.0 and abs(iy - py) <= 1.0, (rot, ix, iy, px, py)
    assert np.isclose(np.nansum(model_in_frame(model, 37.0, center)), model.sum(), rtol=0.02)


def test_roi_mask_and_scaling(rng):
    shape = (61, 61)
    roi = roi_mask(shape, lrmin=5, lrmax=25, roirange=(40, 12, 30.0))
    assert roi.any() and not roi[30, 30]  # centre excluded by lrmin
    with pytest.raises(ValueError):
        roi_mask(shape, lrmin=100)
    model = rng.normal(size=shape) + 3
    ref = 2.5 * model + rng.normal(scale=0.1, size=shape)
    assert abs(scale_model_ls(ref, model, roi) - 2.5) < 0.02
    sigma = np.full(shape, 0.1)
    assert abs(scale_model_ls(ref, model, roi, sigma) - 2.5) < 0.02
    cm = chisq_map(ref, 2.5 * model, sigma, roi)
    assert cm[~roi].sum() == 0 and 0.5 < cm[roi].mean() < 1.5
    red, n, binned = binned_chisq(cm, 3.0)
    assert n > 0 and 0.3 < red < 3.0 and binned.shape == (20, 20)


def test_chisq_scale_and_covariance(rng):
    model = np.linspace(1, 2, 20)
    obs = 3.0 * model + rng.normal(scale=0.05, size=20)
    e = np.full(20, 0.05)
    r = chisq(obs, e, model)
    assert abs(r.alpha - 3.0) < 0.05 and 0.3 < r.reduced < 3.0 and r.ndof == 19
    r2 = chisq(obs, e, model, exclude=[0, 1])
    assert r2.ndof == 17
    psi = np.eye(20)
    rc = chisq_covar(obs, e, model, psi)
    assert abs(rc.alpha - r.alpha) < 1e-3 and abs(rc.chisq - r.chisq) < 1e-3
    rp = chisq(obs, e, model, phot=[4.5], ephot=[0.1], modelphot=[1.5])
    assert rp.ndof == 20
    rm = chisq(obs, e, model, err_model=np.full(20, 0.02))
    assert abs(rm.alpha - 3.0) < 0.1


def test_gb_model_zero_lag_is_one():
    p = np.array([0.2, 0.3, 0.4, 0.8])
    assert np.isclose(gb_model(np.array([0.0]), p, 10.0)[0], 0.5)
    lam = np.linspace(1.15, 2.4, 22)
    lr = (lam[:, None] - lam[None, :]) / np.median(lam)
    c = covariance_matrix(lr, 10.0, 0.2, 0.3, 0.5, 0.4, 0.8)
    assert np.allclose(np.diag(c), 1.0) and np.allclose(c, c.T)
    assert np.all(np.linalg.eigvalsh(c) > -1e-9)
