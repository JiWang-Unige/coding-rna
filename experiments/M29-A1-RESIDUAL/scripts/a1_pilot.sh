#!/bin/bash
# sequential lr pilots (same draws, same seed): 100 steps x batch 8, constant lr after 10% warmup; lr 0 = offset-only baseline
cd $HOME/coding-rna-reorg/experiments/M29-A1-RESIDUAL/scripts
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
CPUS=2,4,6,8,10,12,14,16-24
for LR in 0 3e-5 1e-4 3e-4; do
  taskset -c $CPUS $HOME/envs/generanno/bin/python a1_train.py smoke --out $HOME/coding-rna/outputs/M29-A1-RESIDUAL/pilot_lr$LR --steps 100 --batch 8 --workers 12 --cpus $CPUS --lr-peak $LR --const-lr > $HOME/coding-rna/logs/a1-pilot-lr$LR.out 2> $HOME/coding-rna/logs/a1-pilot-lr$LR.err < /dev/null
done
