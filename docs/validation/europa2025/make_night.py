"""Build Europa_<night>/py for the Python reduction: raw symlinks, dataset .toml from the 2023
template (target coordinates of a mid-sequence frame, date, frame lists after imprep)."""
import csv, glob, os, re, sys, subprocess
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u
night = sys.argv[1]
B = os.path.expanduser("~/charis-val/datasets/2025")
NIGHTS = {"2025-02-09": (153693, 154085, 154086, 154105), "2025-02-10": (154546, 155160, 155161, 155180),
          "2025-02-11": (155712, 156385, 156386, 156405)}
a, b, sa, sb = NIGHTS[night]
T = f"{B}/Europa_{night}"; py = f"{T}/py"
os.makedirs(f"{py}/data/raw", exist_ok=True); os.makedirs(f"{T}/pdi_flats", exist_ok=True)
for f in glob.glob(f"{B}/raw/{night}/pdi_flats/*_cube.fits"):
    d = f"{T}/pdi_flats/{os.path.basename(f)}"; os.path.lexists(d) or os.symlink(f, d)
raw = sorted(glob.glob(f"{B}/raw/{night}/europa/*_cube.fits"))
for f in raw:
    d = f"{py}/data/raw/{os.path.basename(f)}"; os.path.lexists(d) or os.symlink(f, d)
num = lambda f: int(os.path.basename(f)[4:12])
mid = raw[len(raw) // 2]; h = fits.getheader(mid, 0)
ra = Angle(h["RA"], unit=u.hourangle); dec = Angle(h["DEC"], unit=u.deg)
tpl = open(os.path.expanduser("~/charis-val/datasets/Europa_bb/py/Europa_low.toml")).read()
rep = {r'name = "Europa"': f'name = "Europa"', r"ra_deg = .*": f"ra_deg = {ra.deg:.6f}", r"dec_deg = .*": f"dec_deg = {dec.deg:.6f}",
       r"ra_hms = .*": f'ra_hms = "{ra.to_string(unit=u.hourangle, sep=":", precision=2)}"',
       r"dec_dms = .*": f'dec_dms = "{dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True)}"',
       r'date = "20230730"': f'date = "{night.replace("-", "")}"', r"mjd = .*\n": "",
       r'science = .*': 'science = ""', r'sky = .*': 'sky = ""', r'ladder = .*': 'ladder = ""'}
for k, v in rep.items(): tpl = re.sub(k, v, tpl, count=1)
open(f"{py}/Europa_{night}.toml", "w").write(tpl)
# frame lists (sequence indices of the sorted raw list)
dec_ = {r["filename"]: r["decision"] for r in csv.DictReader(open(f"{B}/frame_decisions_laptop_2026-04-28.csv"))}
sci, sky, rej, undecided = [], [], [], []
for i, f in enumerate(raw, 1):
    n = num(f); fn = os.path.basename(f)
    if sa <= n <= sb: sky.append(i)
    elif a <= n <= b:
        d = dec_.get(fn)
        if d == "accepted": sci.append(i)
        elif d == "rejected": rej.append(i)
        else: undecided.append(n); sci.append(i)
def rng(ix):
    out = []; s = p = None
    for i in ix:
        if s is None: s = p = i
        elif i == p + 1: p = i
        else: out.append(f"{s}-{p}" if s != p else f"{s}"); s = p = i
    if s is not None: out.append(f"{s}-{p}" if s != p else f"{s}")
    return ",".join(out)
open(f"{py}/frames_{night}.txt", "w").write(f"science = \"{rng(sci)}\"\nsky = \"{rng(sky)}\"\n# rejected (laptop decisions 2026-04-28): {rng(rej)}\n# undecided, kept: {undecided}\n")
print(f"{night}: {len(raw)} cubes, science {len(sci)} (rejected {len(rej)}, undecided kept {len(undecided)}), sky {len(sky)}; coords {ra.deg:.5f} {dec.deg:.5f} from CRSA00{num(mid)}")
