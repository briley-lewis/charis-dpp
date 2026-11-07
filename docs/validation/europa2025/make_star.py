"""Set up Stars/<star>_<night>/py for a standard star: raw symlinks, flats of the night, .toml from
the HD 236894 template with the star's name and header coordinates."""
import glob, os, re, sys
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u
night, star = sys.argv[1], sys.argv[2]
B = os.path.expanduser("~/charis-val/datasets/2025")
T = f"{B}/Stars/{star}_{night}"; py = f"{T}/py"
os.makedirs(f"{py}/data/raw", exist_ok=True); os.makedirs(f"{T}/pdi_flats", exist_ok=True)
for f in glob.glob(f"{B}/raw/{night}/pdi_flats/*_cube.fits"):
    d = f"{T}/pdi_flats/{os.path.basename(f)}"; os.path.lexists(d) or os.symlink(f, d)
raw = sorted(glob.glob(f"{B}/raw/{night}/{star}/*_cube.fits"))
if not raw: sys.exit(f"no cubes in {B}/raw/{night}/{star}")
for f in raw:
    d = f"{py}/data/raw/{os.path.basename(f)}"; os.path.lexists(d) or os.symlink(f, d)
h = fits.getheader(raw[len(raw) // 2], 0)
ra = Angle(h["RA"], unit=u.hourangle); dec = Angle(h["DEC"], unit=u.deg)
tpl = open(os.path.expanduser("~/charis-val/datasets/HD236894_bb/py/HD236894_low.toml")).read()
rep = {r'name = "HD236894"': f'name = "{star}"', r"ra_deg = .*": f"ra_deg = {ra.deg:.6f}", r"dec_deg = .*": f"dec_deg = {dec.deg:.6f}",
       r"ra_hms = .*": f'ra_hms = "{ra.to_string(unit=u.hourangle, sep=":", precision=2)}"',
       r"dec_dms = .*": f'dec_dms = "{dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True)}"',
       r"plx_mas = .*": "plx_mas = 0.0", r"distance_pc = .*": "distance_pc = 0.0", r"plx_err_mas = .*": "plx_err_mas = 0.0",
       r"pmra_mas_yr = .*": "pmra_mas_yr = 0.0", r"pmdec_mas_yr = .*": "pmdec_mas_yr = 0.0", r'sptype = .*': 'sptype = "-999"',
       r'date = "20230730"': f'date = "{night.replace("-", "")}"', r"mjd = .*\n": "",
       r'science = .*': 'science = ""', r'sky = .*': 'sky = ""', r'ladder = .*': 'ladder = ""'}
for k, v in rep.items(): tpl = re.sub(k, v, tpl, count=1)
tpl = re.sub(r"(\[target.mags\]\n)J = .*\nH = .*\nK = .*", r"\1J = 0.0\nH = 0.0\nK = 0.0", tpl)
open(f"{py}/{star}_{night}.toml", "w").write(tpl)
print(f"{star} {night}: {len(raw)} cubes; coords {ra.deg:.5f} {dec.deg:.5f}; EXPTIME {h.get('EXPTIME')} LYOT {h.get('X_LYOT')}")
