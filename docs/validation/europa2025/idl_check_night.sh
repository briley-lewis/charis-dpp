#!/bin/bash
# New-IDL cross-check of one 2025 night against the Python run: imprep/flat/sky on the same cubes
# (port-mode comparison), then hwp_match and charis_pdi on the Python registered frames with the
# Python (cycle-filtered) hwp_info.txt, Mueller on and off, and the comparisons.
set -u
source ~/anaconda3/etc/profile.d/conda.sh; conda activate charisenv
N=$1; B=~/charis-val/datasets/2025; E=$B/Europa_$N
NEW="+$HOME/charis-val/new-idl/idl:+$HOME/IDLAstro/pro:<IDL_DEFAULT>"
run_idl () { cd $E/$1; CHARIS_DPP_PATH=$HOME/charis-val/new-idl/idl IDL_PATH="$NEW" ~/idl/idl89/bin/idl -quiet -e "$2" > $3.log 2>&1 < /dev/null; echo "$(date +%H:%M) $1 $3: exit $?  $(grep -vE "^% Compiled|READFITS|READCOL|Licens|License:|charis_set_epoch|CHECK_FITS|LAPACK" $3.log | tail -1)"; }
INFO="'Europa_low.info'"
mkdir -p $E/new/data/raw $E/new/pdi_flats $E/new/log $E/new/reduc/prep $E/new/reduc/proc $E/new/reduc/reg; for f in $E/py/data/raw/*_cube.fits; do ln -sf $(readlink -f $f) $E/new/data/raw/; done
python $B/make_info.py $N || exit 1
cp $E/Europa_low.info.master $E/new/Europa_low.info
run_idl new "charis_imprep, $INFO" imprep
cp $E/new/Europa_low.info $E/new/Europa_low.info.after_imprep; cp $E/Europa_low.info.master $E/new/Europa_low.info
run_idl new "charis_pdi_flatfield, $INFO, brightflatdir='../pdi_flats/'" flatfield
run_idl new "charis_pdi_subtract_sky, $INFO, /scalesky" subtract_sky
O=~/charis-val/compare/Europa2025_${N}_new_vs_py_stage1; mkdir -p $O; cd $E; charis compare-trees new/reduc py/reduc --mode port --out $O > $O/stdout.txt 2>&1; echo "--- stage1 new vs py (port): $(tail -2 $O/stdout.txt | head -1)"
mkdir -p $E/new/reduc/reg; find $E/new/reduc/reg -maxdepth 1 -name "*.fits" -delete; cp $E/py/reduc/reg/*.fits $E/new/reduc/reg/; cp $E/py/reduc_left.log $E/py/reduc_right.log $E/new/; cp $E/py/reduc_left.log $E/new/reduc.log
run_idl new "charis_pdi_hwp_match, $INFO, /use_headers" hwp_match
cp $E/new/hwp_info.txt $E/new/hwp_info_idl.txt
diff <(tail -n +2 $E/new/hwp_info_idl.txt | awk '{print $1,$2,$3,$4}') <(tail -n +2 $E/py/hwp_info_unfiltered.txt | awk '{print $1,$2,$3,$4}') > /dev/null && echo "hwp_info IDL vs Python (unfiltered): identical" || echo "hwp_info IDL vs Python (unfiltered): DIFFER ($(diff <(tail -n +2 $E/new/hwp_info_idl.txt | awk '{print $1,$2,$3,$4}') <(tail -n +2 $E/py/hwp_info_unfiltered.txt | awk '{print $1,$2,$3,$4}') | grep -c '^[<>]') lines)"
cp $E/py/hwp_info.txt $E/new/hwp_info.txt
run_idl new "charis_pdi, $INFO, /meanadd, /save_derot_seq" pdi
M=$E/new_nomm; rm -rf $M/reduc/proc; mkdir -p $M/reduc/proc $M/log $M/data; ln -sfn $E/new/reduc/reg $M/reduc/reg; ln -sfn $E/new/data/raw $M/data/raw; cp $E/new/Europa_low.info $E/new/hwp_info.txt $E/new/reduc*.log $M/
run_idl new_nomm "charis_pdi, $INFO, /meanadd, /save_derot_seq, /no_mm_ip_corr, /no_mm_crosstalk_corr" pdi
cd $E; for pair in "new py" "new_nomm py_nomm"; do set -- $pair; O=~/charis-val/compare/Europa2025_${N}_$1_vs_$2; mkdir -p $O; charis compare-trees $1/reduc $2/reduc --mode port --out $O > $O/stdout.txt 2>&1; echo "--- $1 vs $2 (port): $(tail -2 $O/stdout.txt | head -1)"; grep -E "pdi_(q|u|pi|aolp|qphi)(_collapsed)?\.fits" $O/stdout.txt | awk '{printf "  %-30s %-6s %10s %8s %7s\n", $1,$2,$4,$5,$6}'; done
