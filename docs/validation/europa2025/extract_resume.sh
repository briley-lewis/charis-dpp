#!/bin/bash
# Re-run extraction groups until every frame has a cube or is on the failed list (CHARIS-DEP aborts
# the whole call on a frame whose ramp fit divides by zero, e.g. CRSA00154539).
# Usage: extract_resume.sh            -> all groups;   extract_resume.sh <night> <group> -> one group
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
B=~/charis-val/datasets/2025; INI=$B/charis_2025.ini; S=/home/shared/exoserver/Subaru-Europa
FAILED=$B/extract_failed.txt; touch $FAILED
ONLY_N=${1:-}; ONLY_G=${2:-}
run () { # night group first last
  [ -n "$ONLY_N" ] && { [ "$1" = "$ONLY_N" ] && [ "$2" = "$ONLY_G" ] || return; }
  D=$B/raw/$1/$2; mkdir -p $D; cd $D
  for attempt in 1 2 3 4 5 6 7 8; do
    files=""; for n in $(seq $3 $4); do grep -q "CRSA00$n" $FAILED && continue; f=$S/$1/CRSA00$n.fits; [ -f $f ] || f=$(ls $S/$1/*/CRSA00$n.fits 2>/dev/null | head -1); [ -f "$f" ] && [ ! -f $D/CRSA00${n}_cube.fits ] && files="$files $f"; done
    [ -z "$files" ] && { echo "$(date +%H:%M) $1/$2: complete ($(ls $D/*_cube.fits 2>/dev/null | wc -l) cubes)"; return; }
    echo "$(date +%H:%M) $1/$2: attempt $attempt, $(echo $files | wc -w) frames"
    extractcube $files $INI > extract_resume_$attempt.log 2>&1; rc=$?
    if [ $rc -ne 0 ]; then for f in $files; do n=$(basename $f .fits); [ -f $D/${n}_cube.fits ] || { echo "$n $1/$2 $(date +%F_%H:%M)" >> $FAILED; echo "$(date +%H:%M) $1/$2: $n failed, skipping"; break; }; done; fi
  done
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
echo "$(date +%H:%M) resume finished${ONLY_N:+ ($ONLY_N/$ONLY_G)}; failed list:"; cat $FAILED
