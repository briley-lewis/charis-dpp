"""Median left/right beam peak positions (channel 10) per standard, from the raw cubes."""
import glob
import numpy as np
from astropy.io import fits
B = "/home/blewis/charis-val/datasets/2025"
groups = {"HD43384": "2025-02-09", "HD25532": "2025-02-09", "HD40280": "2025-02-09", "HD30675": "2025-02-10", "HD35532": "2025-02-10", "HD283800": "2025-02-11"}
def peak(im, y0, y1, x0, x1):
    sub = im[y0:y1, x0:x1]; iy, ix = np.unravel_index(np.nanargmax(sub), sub.shape); return ix + x0, iy + y0, sub[iy, ix]
for star, night in groups.items():
    files = sorted(glob.glob(f"{B}/raw/{night}/{star}/*_cube.fits"))
    L, R, fl = [], [], []
    for f in files:
        im = fits.getdata(f, 1)[10]
        lx, ly, lm = peak(im, 55, 110, 40, 100); rx, ry, rm = peak(im, 90, 145, 105, 165)
        L.append((lx, ly)); R.append((rx, ry)); fl.append(lm)
    L, R, fl = np.array(L), np.array(R), np.array(fl)
    print(f"{star} {night}: n={len(files)} left {np.median(L,0)} right {np.median(R,0)} sep {np.median(R-L,0)} mid {np.median((L+R)/2,0)} "
          f"| left scatter {L.std(0).round(1)} | ch10 peak median {np.median(fl):.0f} min {fl.min():.0f} max {fl.max():.0f}")
