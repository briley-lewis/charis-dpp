from __future__ import annotations

import numpy as np
import pytest
from astropy.io.fits import Header
from charis_dpp.geometry.apertures import aper_idl, sumaper_im
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.geometry.masks import annulus, loci_zones, sectors, spider_mask
from charis_dpp.image.centroid import center_of_symmetry, cntrd, fit_gaussian2d, gcntrd, xcorr_shift
from charis_dpp.image.filters import (
    fixpix,
    idl_filter_image,
    median_highpass,
    radial_profile,
    radial_profile_2d,
)
from charis_dpp.image.interp import keys_cubic, rescale_radial, rotate_cube, rotate_image, shift_sub
from charis_dpp.image.snr import snr_map
from charis_dpp.image.stats import (
    combine_frames,
    idl_median,
    resistant_mean,
    robust_poly_fit,
    robust_sigma,
)
from charis_dpp.image.wcs import northup, rotation_from_header


def gauss(shape, x0, y0, fwhm=4.7, amp=1.0):
    y, x = np.mgrid[0 : shape[0], 0 : shape[1]]
    s = fwhm / 2.35482
    return amp * np.exp(-0.5 * ((x - x0) ** 2 + (y - y0) ** 2) / s**2)


def test_image_center_odd_even():
    assert image_center((101, 101)) == (50, 50)
    assert image_center((100, 120)) == (60, 50)
    r = radial_grid((101, 101))
    assert r[50, 50] == 0 and r[50, 60] == 10


def test_keys_cubic_reproduces_quadratics():
    y, x = np.mgrid[0:60, 0:60].astype(float)
    im = 2 * x**2 - 3 * x * y + 0.5 * y**2 + x - 7
    xs = np.array([10.3, 20.7, 33.1])
    ys = np.array([12.2, 40.9, 5.6])
    exact = 2 * xs**2 - 3 * xs * ys + 0.5 * ys**2 + xs - 7
    assert np.allclose(keys_cubic(im, xs, ys), exact, atol=1e-9)
    assert np.isnan(keys_cubic(im, np.array([-1.0]), np.array([3.0]))[0])


def test_shift_roundtrip_and_edges():
    g = gauss((81, 81), 40.3, 39.6)
    s = shift_sub(g, 1.25, -0.75)
    back = shift_sub(s, -1.25, 0.75)
    assert np.nanmax(np.abs(back - g)[6:-6, 6:-6]) < 1e-2
    si = shift_sub(g, 2, 0)
    assert np.isnan(si[:, :2]).all() and np.allclose(si[:, 2:], g[:, :-2])
    assert not np.isnan(shift_sub(g, 2, 0, edge="idl")).any()


def test_rotation_direction_matches_rotat():
    g = gauss((101, 101), 60, 50)
    r = rotate_image(g, 90)
    iy, ix = np.unravel_index(np.nanargmax(r), r.shape)
    assert (ix, iy) == (50, 40)  # rotat forward map: x0 + (yi - y0), y0 - (xi - x0)
    back = rotate_image(rotate_image(g, 37.0), -37.0)
    assert np.nanmax(np.abs(back - g)[10:-10, 10:-10]) < 2e-2
    cube = np.stack([g, 2 * g])
    rc = rotate_cube(cube, 90)
    assert np.allclose(rc[1], 2 * rc[0], equal_nan=True)
    g2 = g.copy()
    g2[50, 60] = np.nan
    assert np.isnan(rotate_image(g2, 0.0)[50, 60])


def test_rescale_radial_roundtrip():
    g = gauss((101, 101), 50.3, 49.6)
    back = rescale_radial(rescale_radial(g, 1.3), 1 / 1.3)
    assert np.nanmax(np.abs(back - g)[20:-20, 20:-20]) < 1e-2


def test_centroids_recover_position(rng):
    g = gauss((81, 81), 40.3, 39.6) + rng.normal(0, 0.005, (81, 81))
    assert np.allclose(gcntrd(g, 42, 38, 4.7), (40.3, 39.6), atol=0.03)
    assert np.allclose(cntrd(g, 42, 38, 4.7), (40.3, 39.6), atol=0.05)
    f = fit_gaussian2d(g, 42, 38, 4.7)
    assert abs(f["x"] - 40.3) < 0.02 and abs(f["fwhm_x"] - 4.7) < 0.1
    assert np.isnan(gcntrd(g, 2, 2, 4.7)[0])  # too close to the edge


def test_center_of_symmetry_and_xcorr():
    sym = (
        gauss((101, 101), 50.3, 49.6)
        + gauss((101, 101), 40.3, 49.6)
        + gauss((101, 101), 60.3, 49.6)
    )
    r = radial_grid((101, 101))
    roi = (r > 4) & (r < 30)
    xc, yc = center_of_symmetry(sym, 51, 49, roi)
    assert abs(xc - 50.3) < 0.02 and abs(yc - 49.6) < 0.02
    ref = gauss((101, 101), 50, 50) + gauss((101, 101), 62, 55, amp=0.4)
    moved = shift_sub(ref, 1.3, -2.2)
    x, y = xcorr_shift(ref, np.nan_to_num(moved), roi & np.isfinite(moved), 0, 0)
    assert abs(x - (50 + 1.3)) < 0.02 and abs(y - (50 - 2.2)) < 0.02


def test_apertures():
    flat = np.ones((51, 51))
    flux, eflux, sky, _ = aper_idl(flat, 25, 25, 5.0)
    assert flux == pytest.approx(np.pi * 25, rel=1e-6) and sky == 0
    flux2, _, sky2, _ = aper_idl(flat + 3.0, 25, 25, 5.0, skyrad=(8, 12), setskyval=None)
    assert flux2 == pytest.approx(0.0, abs=1e-9) and sky2 == pytest.approx(4.0)
    s = sumaper_im(flat, 2.0, 20)
    assert s[25, 25] == pytest.approx(np.pi * 4, rel=1e-6) and s[0, 0] == 0


def test_stats(rng):
    x = rng.normal(0, 1, 20000)
    m, s, nrej = resistant_mean(x, 3.0, return_all=True)
    assert abs(m) < 0.03 and 0.005 < s < 0.012
    assert abs(robust_sigma(x) - 1) < 0.03
    x2 = x.copy()
    x2[:50] = 100.0
    assert abs(resistant_mean(x2, 3.0)) < 0.05
    assert idl_median([1, 2, 3, 4], even=False) == 3 and idl_median([1, 2, 3, 4]) == 2.5
    stack = np.stack([np.full((4, 4), v) for v in (1.0, 2.0, 3.0, 100.0)])
    assert combine_frames(stack, "median")[0, 0] == 2.5
    assert abs(combine_frames(stack, "mean")[0, 0] - 2.0) < 1e-6
    xx = np.linspace(0, 10, 100)
    yy = 1 + 2 * xx + 0.1 * xx**2
    yy[::10] += 50
    c = robust_poly_fit(xx, yy, 2)
    assert np.allclose(c, [1, 2, 0.1], atol=0.05)


def test_filters_and_profiles():
    y, x = np.mgrid[0:101, 0:101]
    r = np.hypot(x - 50, y - 50).astype(float)
    radii, prof = radial_profile(r, 1.0)
    assert radii[1] == 1 and abs(prof[10] - 10) < 0.3
    p2 = radial_profile_2d(r, 1.0, 0, 40)
    assert np.nanmax(np.abs(p2 - r)[r < 39]) < 0.3
    assert p2[50, 95] == 0  # outside rmax untouched
    im = np.zeros((41, 41))
    im[20, 20] = 10
    f = idl_filter_image(im, median=5)
    assert f[20, 20] == 0 and np.array_equal(f[:2, :], im[:2, :])
    hp = median_highpass(im, 5)
    assert hp[20, 20] == 10


def test_fixpix_single_hot_pixel(rng):
    im = rng.normal(0, 1, (60, 60))
    im[30, 31] += 200.0
    im[14, 22] -= 200.0
    fixed, mask = fixpix(im, 3, return_mask=True)
    assert mask[30, 31] == 1 and mask[14, 22] == -1
    assert abs(fixed[30, 31]) < 5 and abs(fixed[14, 22]) < 5


def test_masks_and_zones():
    a = annulus((101, 101), 10, 15)
    r = radial_grid((101, 101))
    assert a.sum() == ((r >= 10) & (r < 15)).sum()
    s = sectors((101, 101), 4)
    assert set(np.unique(s)) == {0, 1, 2, 3}
    zones = loci_zones((101, 101), 10, 30, 5, 3.0, 50, 1.0)
    assert len(zones) > 4
    total = sum(z.sub.sum() for z in zones)
    assert total == ((r >= 10) & (r < 30)).sum()  # zones tile the annuli
    for z in zones:
        assert z.opt.sum() >= z.sub.sum() * 0.5
    sp = spider_mask((101, 101), 3.0, [0, 90, 180, 270])
    assert sp[50, 80] == 0 and sp[50, 50] == 0 and sp[80, 80] == 1
    assert sp[50 + 2, 80] > 0  # fractional edge
    assert np.isclose(sp[53, 80], 1.0)


def test_snr_map_detects_source(rng):
    im = rng.normal(0, 1, (101, 101))
    im += gauss((101, 101), 71, 50, fwhm=3.0, amp=40)
    res = snr_map(im, 3.0, 45, coord=(71, 50))
    assert res.snr > 20
    assert abs(res.position[0] - 71) < 0.3
    res2 = snr_map(im, 3.0, 45, finite=True)
    assert np.nanmax(res2.snrmap[45:55, 65:77]) > 10


def test_northup_updates_wcs_once():
    g = gauss((101, 101), 60, 50)
    h = Header()
    a = np.deg2rad(30.0)
    s = 0.0162 / 3600
    h["CD1_1"], h["CD1_2"], h["CD2_1"], h["CD2_2"] = (
        np.cos(a) * s,
        np.sin(a) * s,
        -np.sin(a) * s,
        np.cos(a) * s,
    )
    h["CRPIX1"], h["CRPIX2"] = 51, 51
    assert abs(rotation_from_header(h) - 30.0) < 1e-9
    out, rotang = northup(g, h, 2.0)
    assert rotang == pytest.approx(28.0)
    assert h["ROTANG"] == pytest.approx(28.0)
    assert abs(rotation_from_header(h) - 2.0) < 1e-9  # only the north offset remains
    iy, ix = np.unravel_index(np.nanargmax(out), out.shape)
    t = np.deg2rad(28.0)
    assert abs(ix - (50 + 10 * np.cos(t))) < 1.1 and abs(iy - (50 - 10 * np.sin(t))) < 1.1
