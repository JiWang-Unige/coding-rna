#!/bin/bash
# Stage 1 (sequential, resumable): prep -> reference check -> alignment -> GPU inference of B.   Marker: $O/stage1.done
O=$HOME/coding-rna/outputs/M37-JOINT-DECODE-SMOKE; mkdir -p $O; L=$HOME/coding-rna/logs
S=$HOME/coding-rna-reorg/experiments/M37-JOINT-DECODE-SMOKE/scripts
PY=$HOME/envs/generanno/bin/python
cd $S
[ -f $HOME/coding-rna/data/m37/oryza_sativa_9311/source.json ] || nice -n 10 taskset -c 25-31 $PY m37_prep.py > $L/m37-prep.out 2> $L/m37-prep.err || { echo prep_failed >> $O/stage1.log; exit 1; }
echo "prep ok $(date -Is)" >> $O/stage1.log
[ -f $O/reference_counts.json ] || nice -n 10 taskset -c 25-31 $PY m37_check_ref.py > $L/m37-checkref.out 2> $L/m37-checkref.err
echo "checkref rc=$? $(date -Is)" >> $O/stage1.log
[ -s $O/align/A_vs_B.paf ] || bash m37_align.sh
echo "align done $(date -Is)" >> $O/stage1.log
ls $O/inferB/DONE_*.json >/dev/null 2>&1 || bash m37_gpu_infer.sh
echo "infer rc=$? $(date -Is)" >> $O/stage1.log
touch $O/stage1.done
