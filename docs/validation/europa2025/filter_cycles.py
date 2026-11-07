"""Drop HWP cycles whose four frames span more than `maxspan` minutes (culled frames can leave
stragglers that the nearest-in-time matcher pairs across gaps). Rewrites hwp_info.txt in place."""
import re, sys, glob, os, numpy as np
from astropy.io import fits
from charis_dpp.pdi.hwp_match import HwpTable
path = sys.argv[1]; maxspan = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
t = HwpTable.load(path)
mjd = {}
for f in glob.glob("reduc/prep/n*e.fits"):
    h = fits.getheader(f, 0); n = int(os.path.basename(f)[1:5]); mjd[n] = float(h.get("MJD", fits.getheader(f, 1).get("MJD", np.nan)))
names = np.array([int(re.search(r"\d+", str(x)).group()) for x in t.names])  # "12" or "n0012leftreg.fits"
tm = np.array([mjd.get(n, np.nan) for n in names]) * 1440.0
spans = {}
for c in sorted(set(t.hwp_cycle[t.hwp_cycle >= 0])):
    m = t.hwp_cycle == c; spans[c] = float(np.nanmax(tm[m]) - np.nanmin(tm[m]))
bad = [c for c, s in spans.items() if s > maxspan]
for c in bad: t.hwp_cycle[t.hwp_cycle == c] = -1
# renumber the surviving cycles consecutively
keep = sorted(set(t.hwp_cycle[t.hwp_cycle >= 0])); remap = {c: i for i, c in enumerate(keep)}
t.hwp_cycle = np.array([remap.get(c, -1) for c in t.hwp_cycle])
t.save(path)
print(f"{path}: {len(spans)} cycles, median span {np.median(list(spans.values())):.2f} min, dropped {len(bad)} with span > {maxspan} min: " + ", ".join(f"{c}({spans[c]:.1f})" for c in bad))
