#!/bin/bash
# Carbon-A native post-processing (CPU numba Viterbi + confidence network) on saved raw probabilities -> prediction.gff (all genes with confidence=)
# usage: m30_native.sh <species> <contig>
SP=$1; C=$2
D=$HOME/coding-rna/outputs/M30-HOLDOUT-ADAPT/carbonA_raw/$SP/$C
cd $HOME/coding-rna-reorg/experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/src
export PYTHONPATH=$HOME/envs/m29_extra OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
[ -f $D/prediction.gff ] && exit 0
taskset -c 18-23 $HOME/envs/generanno/bin/python post_processing.py --input $(ls $D/*.raw_probs.parquet) --output $D/prediction --viterbi_backend numba --device cpu > $HOME/coding-rna/logs/m30-native-$C.out 2> $HOME/coding-rna/logs/m30-native-$C.err < /dev/null
ls $D
