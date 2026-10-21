#!/bin/bash
# usage: c1_infer_retry.sh <outdir> <model.pt> <regions> [decode_workers]   retrying, resumable per (region, strand)
OUT=$1; MODEL=$2; REGIONS=$3; DW=${4:-4}
CPUS=${M32_CPUS:-2,4,6,8,10,12,14,16-24}
PY=$HOME/envs/plantcad2/bin/python
cd $HOME/coding-rna-reorg/experiments/M32-C1-PLANTCAD2/scripts
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
mkdir -p $OUT
for a in $(seq 1 30); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt 5000 ]; then sleep 120; continue; fi
  taskset -c $CPUS $PY c1_infer.py --model $MODEL --regions $REGIONS --out $OUT --decode-workers $DW --cpus $CPUS >> $HOME/coding-rna/logs/c1-infer-$(basename $OUT).out 2>> $HOME/coding-rna/logs/c1-infer-$(basename $OUT).err < /dev/null
  rc=$?
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then echo "{\"event\":\"done\",\"time\":\"$(date -Is)\"}" >> $OUT/interruptions.jsonl; break; fi
  echo "{\"event\":\"interruption\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\",\"stderr_tail\":\"$(tail -n 2 $HOME/coding-rna/logs/c1-infer-$(basename $OUT).err | tr '\n"' '  ' | cut -c1-300)\"}" >> $OUT/interruptions.jsonl
  sleep 30
done
