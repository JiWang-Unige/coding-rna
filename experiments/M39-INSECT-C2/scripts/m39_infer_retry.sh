#!/bin/bash
# usage: m39_infer_retry.sh <outdir> <model.pt> <regions> [decode_workers]   retrying, resumable per (region, strand); waits while the GPU has < 5 GB free
OUT=$1; MODEL=$2; REGIONS=$3; DW=${4:-4}
CPUS=${M39_CPUS:-2,4,6,10,12,14,16-24}   # CPU 8 dropped 2026-11-12 (segfaults cluster on CPU 0 and 8, P6 dmesg)
cd $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts
source ./m39_env.sh
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
mkdir -p $OUT $HOME/coding-rna/logs
NAME=m39-infer-$(basename $(dirname $OUT))-$(basename $OUT)
for a in $(seq 1 40); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt 5000 ]; then sleep 120; continue; fi
  taskset -c $CPUS $PY m39_infer.py --model $MODEL --regions $REGIONS --out $OUT --decode-workers $DW --cpus $CPUS >> $HOME/coding-rna/logs/$NAME.out 2>> $HOME/coding-rna/logs/$NAME.err < /dev/null
  rc=$?
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then echo "{\"event\":\"done\",\"time\":\"$(date -Is)\"}" >> $OUT/interruptions.jsonl; exit 0; fi
  echo "{\"event\":\"interruption\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\",\"stderr_tail\":\"$(tail -n 2 $HOME/coding-rna/logs/$NAME.err | tr '\n"' '  ' | cut -c1-300)\"}" >> $OUT/interruptions.jsonl
  sleep 30
done
exit 9
