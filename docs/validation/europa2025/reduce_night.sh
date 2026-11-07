#!/bin/bash
# Python (charis-dpp) PDI reduction of one Feb 2025 Europa night: imprep, lamp flats, scaled sky
# subtraction, xcorr registration, header HWP matching, charis pdi with and without the Mueller
# correction. Usage: reduce_night.sh 2025-02-10
set -u
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
N=$1; B=~/charis-val/datasets/2025; T=$B/Europa_$N; C=Europa_$N.toml
python $B/make_night.py $N || exit 1
cd $T/py
echo "$(date +%H:%M) imprep"; charis imprep $C --pdi > log_imprep.txt 2>&1 || { tail -3 log_imprep.txt; exit 1; }
python $B/set_frames.py $C frames_$N.txt
echo "$(date +%H:%M) flatfield"; charis pdi-flatfield $C --brightflatdir ../pdi_flats/ > log_flatfield.txt 2>&1 || { tail -3 log_flatfield.txt; exit 1; }
echo "$(date +%H:%M) subtract-sky"; charis pdi-subtract-sky $C --scalesky > log_skysub.txt 2>&1 || { tail -3 log_skysub.txt; exit 1; }
echo "$(date +%H:%M) register (xcorr)"; charis pdi-register $C --method xcorr > log_register.txt 2>&1 || { tail -3 log_register.txt; exit 1; }
echo "$(date +%H:%M) hwp-match"; charis pdi-hwp-match $C --source header > log_hwp.txt 2>&1 || { tail -3 log_hwp.txt; exit 1; }
cp hwp_info.txt hwp_info_unfiltered.txt; python $B/filter_cycles.py hwp_info.txt 10
echo "$(date +%H:%M) pdi (Mueller on)"; charis pdi $C --meanadd --save-derot-seq > log_pdi.txt 2>&1 || { tail -3 log_pdi.txt; exit 1; }
M=$T/py_nomm; mkdir -p $M/reduc/reg $M/reduc/proc $M/data; ln -sfn $T/py/data/raw $M/data/raw; cp $C hwp_info.txt $M/
for f in $T/py/reduc/reg/*.fits; do ln -sf $f $M/reduc/reg/; done
cd $M; echo "$(date +%H:%M) pdi (Mueller off)"; charis pdi $C --meanadd --save-derot-seq --no-mm-ip-corr --no-mm-crosstalk-corr > log_pdi.txt 2>&1 || { tail -3 log_pdi.txt; exit 1; }
echo "$(date +%H:%M) $N done: $(ls $T/py/reduc/proc/pdi_*_collapsed.fits | wc -l) collapsed products"
