"""After imprep (which rewrites the frame lists), put the accepted science and the sky indices back."""
import re, sys
toml, frames = sys.argv[1], sys.argv[2]
s = open(toml).read()
for line in open(frames):
    if line.startswith("science") or line.startswith("sky"):
        k, v = line.split(" = ", 1); s = re.sub(rf'^{k} = .*$', f'{k} = {v.strip()}', s, count=1, flags=re.M)
open(toml, "w").write(s); print(open(frames).read())
