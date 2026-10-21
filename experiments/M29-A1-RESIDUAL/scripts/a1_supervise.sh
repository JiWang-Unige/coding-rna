#!/bin/bash
# Crash-resilient launcher for a1_train.py: restarts from the latest checkpoint and records every interruption.
# usage: a1_supervise.sh <name> <outdir> <steps> <batch> <workers>      (run with nohup; prints nothing)
NAME=$1; OUT=$2; STEPS=$3; BATCH=$4; WORKERS=$5
CPUS=2,4,6,8,10,12,14,16-24
TS=$(date +%Y%m%d-%H%M%S)
LOG=$HOME/coding-rna/logs/$NAME-$TS
mkdir -p "$OUT"
cd $HOME/coding-rna-reorg/experiments/M29-A1-RESIDUAL/scripts
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
for attempt in $(seq 1 60); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt 9000 ]; then
    echo "{\"event\":\"gpu_busy_wait\",\"free_mib\":$free,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; sleep 120; continue
  fi
  taskset -c $CPUS $HOME/envs/generanno/bin/python a1_train.py run --out "$OUT" --steps $STEPS --batch $BATCH --workers $WORKERS --cpus $CPUS --resume >> $LOG.out 2>> $LOG.err < /dev/null &
  PID=$!
  echo "{\"event\":\"launch\",\"attempt\":$attempt,\"pid\":$PID,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  wait $PID; rc=$?
  if [ -f "$OUT/DONE.json" ]; then echo "{\"event\":\"done\",\"rc\":$rc,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; break; fi
  last=$(~/envs/generanno/bin/python -c "import torch;print(torch.load('$OUT/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo none)
  err=$(tail -n 3 $LOG.err | tr '\n"' '  ' | cut -c1-400)
  seg=$(dmesg 2>/dev/null | grep -c segfault)
  echo "{\"event\":\"interruption\",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"dmesg_segfaults_total\":$seg,\"stderr_tail\":\"$err\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  pkill -P $$ -f a1_dp 2>/dev/null; sleep 30
done
