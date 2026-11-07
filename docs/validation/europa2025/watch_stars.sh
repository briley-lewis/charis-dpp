#!/bin/bash
# After the whole queue and the gap-filling pass: reduce the polarised and unpolarised standards.
B=~/charis-val/datasets/2025; cd $B
until grep -q "resume finished" extract_resume.out 2>/dev/null; do sleep 300; done
for ns in "2025-02-09 HD43384" "2025-02-10 HD30675" "2025-02-11 HD283800" "2025-02-09 HD25532" "2025-02-09 HD40280" "2025-02-10 HD35532"; do set -- $ns; ./reduce_star.sh $1 $2 > reduce_star_$2_$1.out 2>&1; tail -3 reduce_star_$2_$1.out; done
echo "$(date +%H:%M) watch_stars finished"
