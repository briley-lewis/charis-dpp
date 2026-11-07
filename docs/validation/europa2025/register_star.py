"""PSF-centroid PDI registration of a standard star with per-beam guesses measured from the data.

`charis pdi-register --method psf` aims gcntrd's 9 px box at frame centre +/- the nominal Wollaston
offset (31.8, 16.5); in Feb 2025 the stars sit up to 7 px from that nominal spot (pointing), so most
frames of HD 43384 / HD 40280 found no centroid. This measures the median channel-10 peak of each beam
over the prepared cubes and offsets the guesses accordingly (the target stays the frame centre).
Usage: python register_star.py <config.toml>
"""
import glob, sys
import numpy as np
from astropy.io import fits
import charis_dpp.pdi.register as R
from charis_dpp.config import DatasetConfig

cfg = DatasetConfig.load(sys.argv[1])
files = sorted(glob.glob("reduc/prep/n*e_flat.fits"))
def peak(im, y0, y1, x0, x1):
    sub = im[y0:y1, x0:x1]; iy, ix = np.unravel_index(np.nanargmax(sub), sub.shape); return ix + x0, iy + y0
L, Rr = [], []
for f in files:
    im = fits.getdata(f, 1)[10]
    L.append(peak(im, 55, 110, 40, 100)); Rr.append(peak(im, 90, 145, 105, 165))
L, Rr = np.median(L, 0), np.median(Rr, 0)
ny, nx = im.shape; c = ((nx - 1) / 2.0, (ny - 1) / 2.0)
goff = {"leftreg": (L[0] - c[0], L[1] - c[1]), "rightreg": (Rr[0] - c[0], Rr[1] - c[1])}
print(f"beam peaks (ch10, median of {len(files)} frames): left {L} right {Rr}; guess offsets {goff}")
orig = R.register_cube
def patched(cfg, **kw):
    kw["guessoffsets"] = goff[kw["out_suffix"]]
    return orig(cfg, **kw)
R.register_cube = patched
rl_l, rl_r = R.pdi_register_cube(cfg, method="psf")
print(f"registered: {len(glob.glob('reduc/reg/n*leftreg.fits'))} left, {len(glob.glob('reduc/reg/n*rightreg.fits'))} right of {len(files)} frames")
