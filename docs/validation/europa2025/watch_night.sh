#!/bin/bash
# Wait for a night's Europa cubes, fill any gaps, then run the Python reduction, stage 3 and the IDL check.
N=$1; B=~/charis-val/datasets/2025; cd $B
until grep -q "$N/europa: done" extract_queue.out 2>/dev/null; do sleep 120; done
./extract_resume.sh $N pdi_flats; ./extract_resume.sh $N europa
./reduce_night.sh $N > reduce_$N.out 2>&1; tail -2 reduce_$N.out
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv; python europa2025_stage3.py $N > stage3_$N.out 2>&1
./idl_check_night.sh $N > idlcheck_$N.out 2>&1
echo "$(date +%H:%M) watch_night $N finished"
