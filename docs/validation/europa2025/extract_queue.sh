#!/bin/bash
# CHARIS-DEP extraction queue for the Feb 2025 nights, one calibration (2025-02-09, CRSA00153676,
# no darks) as on the laptop; cubes land in ~/charis-val/datasets/2025/raw/<night>/<group>/.
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
INI=~/charis-val/datasets/2025/charis_2025.ini
S=/home/shared/exoserver/Subaru-Europa
run () { # night group first last
  D=~/charis-val/datasets/2025/raw/$1/$2; mkdir -p $D; cd $D
  files=""; for n in $(seq $3 $4); do f=$S/$1/CRSA00$n.fits; [ -f $f ] || f=$(ls $S/$1/*/CRSA00$n.fits 2>/dev/null | head -1); [ -f "$f" ] && [ ! -f $D/CRSA00${n}_cube.fits ] && files="$files $f"; done
  [ -z "$files" ] && { echo "$(date +%H:%M) $1/$2: nothing to do"; return; }
  echo "$(date +%H:%M) $1/$2: extracting $(echo $files | wc -w) frames"
  extractcube $files $INI > extract.log 2>&1; rc=$?
  echo "$(date +%H:%M) $1/$2: done, $(ls $D/*_cube.fits 2>/dev/null | wc -l) cubes, exit $rc"
}
run 2025-02-10 pdi_flats 154537 154541
run 2025-02-10 europa    154546 155180
run 2025-02-11 pdi_flats 155698 155702
run 2025-02-11 europa    155712 156405
run 2025-02-09 pdi_flats 153678 153682
run 2025-02-09 europa    153693 154105
run 2025-02-09 HD43384   154405 154484
run 2025-02-10 HD30675   155277 155377
run 2025-02-11 HD283800  156406 156446
run 2025-02-09 HD25532   154118 154218
run 2025-02-09 HD40280   154220 154324
run 2025-02-10 HD35532   155480 155580
echo "$(date +%H:%M) queue finished"
