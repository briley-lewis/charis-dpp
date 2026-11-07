#!/bin/bash
# Reduce the six Feb 2025 standards in parallel (after the registration fixes), then summarise.
B=~/charis-val/datasets/2025; cd $B
for ns in "2025-02-09 HD43384" "2025-02-10 HD30675" "2025-02-11 HD283800" "2025-02-09 HD25532" "2025-02-09 HD40280" "2025-02-10 HD35532"; do set -- $ns; ./reduce_star.sh $1 $2 > reduce_star_$2_$1.out 2>&1 & done
wait
for f in reduce_star_*.out; do echo "== $f"; tail -12 $f; done
echo "$(date +%H:%M) run_stars finished"
