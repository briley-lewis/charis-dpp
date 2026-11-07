#!/bin/bash
# Python PDI reduction of a standard star (no sky frames, PSF-centroid registration with per-beam
# guesses from register_star.py), Mueller on and off. Usage: reduce_star.sh <night> <star>
set -u
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
N=$1; ST=$2; B=~/charis-val/datasets/2025; T=$B/Stars/${ST}_$N; C=${ST}_$N.toml
rm -rf $T/py/reduc $T/py_nomm $T/py/*.log $T/py/log_*.txt $T/py/hwp_info*.txt
python $B/make_star.py $N $ST || exit 1
cd $T/py
echo "$(date +%H:%M) $ST imprep"; charis imprep $C --pdi > log_imprep.txt 2>&1 || { tail -3 log_imprep.txt; exit 1; }
grep -E "^(science|sky) =" $C
echo "$(date +%H:%M) $ST flatfield"; charis pdi-flatfield $C --brightflatdir ../pdi_flats/ > log_flatfield.txt 2>&1 || { tail -3 log_flatfield.txt; exit 1; }
echo "$(date +%H:%M) $ST register (psf, measured beam guesses)"; python $B/register_star.py $C > log_register.txt 2>&1 || { tail -3 log_register.txt; exit 1; }
grep -E "beam peaks|registered:" log_register.txt; grep -c "registration failed" log_register.txt | sed "s/^/frames skipped: /"
echo "$(date +%H:%M) $ST hwp-match"; charis pdi-hwp-match $C --source header > log_hwp.txt 2>&1 || { tail -3 log_hwp.txt; exit 1; }
cp hwp_info.txt hwp_info_unfiltered.txt; python $B/filter_cycles.py hwp_info.txt 10
echo "$(date +%H:%M) $ST pdi (Mueller on)"; charis pdi $C --meanadd --save-derot-seq > log_pdi.txt 2>&1 || { tail -3 log_pdi.txt; exit 1; }
M=$T/py_nomm; mkdir -p $M/reduc/reg $M/reduc/proc $M/data; ln -sfn $T/py/data/raw $M/data/raw; cp $C hwp_info.txt $M/
for f in $T/py/reduc/reg/*.fits; do ln -sf $f $M/reduc/reg/; done
cd $M; echo "$(date +%H:%M) $ST pdi (Mueller off)"; charis pdi $C --meanadd --save-derot-seq --no-mm-ip-corr --no-mm-crosstalk-corr > log_pdi.txt 2>&1 || { tail -3 log_pdi.txt; exit 1; }
cd $T; python $B/star_pol.py py py_nomm
echo "$(date +%H:%M) $ST finished"
