#!/bin/bash
# inference-only diagnostic on the INTERNAL validation chromosomes: shrink the saved head logits by alpha (re-decode, CPU only)
O=$HOME/coding-rna/outputs/M29-A2-TEMPERATURE
cd $HOME/coding-rna-reorg/experiments/M29-A2-TEMPERATURE/scripts
export PYTHONUNBUFFERED=1
for TAU in 0.05 0.03; do for A in ${ALPHAS:-0.25 0.5}; do
  OUT=$O/valshrink_tau${TAU}_a$A
  taskset -c 2,4,6,8,10,12,14,16-24 $HOME/envs/generanno/bin/python $HOME/coding-rna-reorg/experiments/M29-A1-RESIDUAL/scripts/a1_infer.py --head-from $O/val_tau${TAU}_256 --head-scale $A --seqs NC_003071.7,NC_089044.1 --out $OUT --decode-workers 4 --cpus 2,4,6,8,10,12,14,16-24 > $HOME/coding-rna/logs/a2-shrink-$TAU-$A.out 2>&1 < /dev/null
  taskset -c 18-23 $HOME/envs/generanno/bin/python a2_valscore.py $OUT 2>&1 | grep -v "Warn\|BAOBAB" | cut -c1-500 >> $O/shrink_diag.log
done; done
