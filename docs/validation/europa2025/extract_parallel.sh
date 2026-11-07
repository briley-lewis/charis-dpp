#!/bin/bash
# Parallel CHARIS-DEP extraction of the remaining Feb 2025 groups (extractcube is effectively
# single-threaded, ~22 s/frame, 3.7 GB). Each chunk reuses the extract_resume.sh retry logic on its own
# frame range with the shared failed list. Signals: "<night>/europa: done" -> extract_queue.out
# (watch_night.sh), "resume finished" -> extract_resume.out (watch_stars.sh).
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
B=~/charis-val/datasets/2025; INI=$B/charis_2025.ini; S=/home/shared/exoserver/Subaru-Europa
FAILED=$B/extract_failed.txt; touch $FAILED
RUNNING_0210=${1:-}
chunk () { # night group first last
  D=$B/raw/$1/$2; mkdir -p $D; cd $D
  for attempt in 1 2 3 4 5 6 7 8; do
    files=""; for n in $(seq $3 $4); do grep -q "CRSA00$n" $FAILED && continue; f=$S/$1/CRSA00$n.fits; [ -f $f ] || f=$(ls $S/$1/*/CRSA00$n.fits 2>/dev/null | head -1); [ -f "$f" ] && [ ! -f $D/CRSA00${n}_cube.fits ] && files="$files $f"; done
    [ -z "$files" ] && { echo "$(date +%H:%M) $1/$2 $3-$4: complete"; return; }
    echo "$(date +%H:%M) $1/$2 $3-$4: attempt $attempt, $(echo $files | wc -w) frames"
    extractcube $files $INI > extract_$3-$4_$attempt.log 2>&1; rc=$?
    if [ $rc -ne 0 ]; then for f in $files; do n=$(basename $f .fits); [ -f $D/${n}_cube.fits ] || { echo "$n $1/$2 $(date +%F_%H:%M)" >> $FAILED; echo "$(date +%H:%M) $1/$2 $3-$4: $n failed, skipping"; break; }; done; fi
  done
}
# the sequential queue's running 02-10 job: signal its watcher when it exits
if [ -n "$RUNNING_0210" ]; then
  ( while kill -0 $RUNNING_0210 2>/dev/null; do sleep 60; done
    echo "$(date +%H:%M) 2025-02-10/europa: done, $(ls $B/raw/2025-02-10/europa/*_cube.fits | wc -l) cubes (queue job finished; parallel mode)" >> $B/extract_queue.out ) &
fi
P11=(); P09=(); PS=()
chunk 2025-02-11 pdi_flats 155698 155702 & P11+=($!)
for r in 155712-155850 155851-155989 155990-156128 156129-156267 156268-156405; do chunk 2025-02-11 europa ${r%-*} ${r#*-} & P11+=($!); done
chunk 2025-02-09 pdi_flats 153678 153682 & P09+=($!)
for r in 153693-153830 153831-153968 153969-154105; do chunk 2025-02-09 europa ${r%-*} ${r#*-} & P09+=($!); done
chunk 2025-02-09 HD43384  154405 154484 & PS+=($!)
chunk 2025-02-10 HD30675  155277 155377 & PS+=($!)
chunk 2025-02-11 HD283800 156406 156446 & PS+=($!)
chunk 2025-02-09 HD25532  154118 154218 & PS+=($!)
chunk 2025-02-09 HD40280  154220 154324 & PS+=($!)
chunk 2025-02-10 HD35532  155480 155580 & PS+=($!)
echo "$(date +%H:%M) launched $(( ${#P11[@]} + ${#P09[@]} + ${#PS[@]} )) chunks"
wait ${P11[@]}; echo "$(date +%H:%M) 2025-02-11/europa: done, $(ls $B/raw/2025-02-11/europa/*_cube.fits | wc -l) cubes (parallel)" | tee -a $B/extract_queue.out
wait ${P09[@]}; echo "$(date +%H:%M) 2025-02-09/europa: done, $(ls $B/raw/2025-02-09/europa/*_cube.fits | wc -l) cubes (parallel)" | tee -a $B/extract_queue.out
wait ${PS[@]};  echo "$(date +%H:%M) stars done; resume finished (parallel); failed list:" | tee -a $B/extract_resume.out; cat $FAILED
echo "$(date +%H:%M) extract_parallel finished"
