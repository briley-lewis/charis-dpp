"""IDL .info for a 2025 Europa night from the 2023 template (docs/validation/Europa_low.info):
same cards, with the coordinates/date of the Python run and the Python frame lists."""
import os, re, sys, warnings
from astropy.io import fits
warnings.simplefilter("ignore")
night = sys.argv[1]; B = os.path.expanduser("~/charis-val/datasets/2025"); T = f"{B}/Europa_{night}"
toml = open(f"{T}/py/Europa_{night}.toml").read()
g = lambda k: re.search(rf"^{k} = (.*)$", toml, re.M).group(1).strip().strip('"')
tpl = open(os.path.expanduser("~/charis-val/new-idl/docs/validation/Europa_low.info")).read().splitlines()
vals = {"NAME": "Europa", "RA": float(g("ra_deg")), "DEC": float(g("dec_deg")), "RA_HMS": g("ra_hms"), "DEC_DMS": g("dec_dms"),
        "OBSDATE": night.replace("-", ""), "FNUM_SAT": g("science"), "FNUM_SKY": g("sky")}
out = []
for line in tpl:
    key = line[:8].strip()
    if key in vals:
        c = fits.Card.fromstring(line); comment = c.comment; c.value = vals[key]; line = c.image
        if isinstance(vals[key], str) and len(vals[key]) > 68:
            # the IDL `param`/sxpar_param reader takes the value from the whole (unlimited) line but
            # does not follow CONTINUE cards without LONGSTRN: write one long line, as a hand edit would
            line = f"{key:<8}= '{vals[key]}' / {comment}"
    out.append(line)
open(f"{T}/Europa_low.info.master", "w").write("\n".join(out) + "\n")
print("\n".join(l for l in out if l[:8].strip() in vals))
