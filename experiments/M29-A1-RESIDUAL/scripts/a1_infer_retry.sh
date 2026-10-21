#!/bin/bash
# retrying wrapper for a1_infer.py (resumable per contig/strand); logs every interruption
OUT=$1; MODEL=$2; SEQS=$3
CPUS=2,4,6,8,10,12,14,16-24
cd $HOME/coding-rna-reorg/experiments/M29-A1-RESIDUAL/scripts
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
mkdir -p $OUT
for a in $(seq 1 20); do
  taskset -c $CPUS $HOME/envs/generanno/bin/python a1_infer.py --model $MODEL --seqs $SEQS --out $OUT --decode-workers 6 --cpus $CPUS >> $HOME/coding-rna/logs/a1-infer-$(basename $OUT).out 2>> $HOME/coding-rna/logs/a1-infer-$(basename $OUT).err < /dev/null
  rc=$?
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then echo "{\"event\":\"done\",\"time\":\"$(date -Is)\"}" >> $OUT/interruptions.jsonl; break; fi
  echo "{\"event\":\"interruption\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\",\"stderr_tail\":\"$(tail -n 2 $HOME/coding-rna/logs/a1-infer-$(basename $OUT).err | tr '\n"' '  ' | cut -c1-300)\"}" >> $OUT/interruptions.jsonl
  sleep 30
done
