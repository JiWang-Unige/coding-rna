#!/bin/bash
# usage: m39_supervise.sh <arm> <outdir> <stop_at> [extra m39_train_lora args...]   (copy of m34_supervise.sh)
# crash-resilient trainer launcher (as M31 c1_supervise_lora.sh): resumes from ckpt_latest.pt after any interruption (GPU OOM from other users, segfault, ...), logs each one in
# interruptions.jsonl, waits while the GPU has < 7.5 GB free, and gives up (FATAL) when 4 consecutive attempts make no progress (a deterministic code error).
ARM=$1; OUT=$2; STOP=$3; shift 3
CPUS=${M39_CPUS:-2,4,6,10,12,14,16-24}   # CPU 8 dropped 2026-11-12 (segfaults cluster on CPU 0 and 8, P6 dmesg)
cd $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts
source ./m39_env.sh
LOG=$HOME/coding-rna/logs/m39-train-$ARMPFX$ARM-$(date +%Y%m%d-%H%M%S)
mkdir -p "$OUT" $HOME/coding-rna/logs
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
laststep=-1; noprog=0
for attempt in $(seq 1 80); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  MINFREE=7500; [ "$M39_BACKEND" = pc2 ] && MINFREE=6000; [ "$M39_BACKEND" = onehot ] && MINFREE=3000          # pc2 arms use ~5 GB (as m34_supervise.sh since 11-07); NT-v2 unchanged
  if [ "$free" -lt $MINFREE ]; then echo "{\"event\":\"gpu_busy_wait\",\"free_mib\":$free,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; sleep 120; continue; fi
  taskset -c $CPUS $PY m39_train_lora.py run --arm $ARM --out "$OUT" --steps 3072 --stop-at $STOP --batch 8 --workers 12 --cpus $CPUS --resume ${M39_LORA:+--lora $M39_LORA} "$@" >> $LOG.out 2>> $LOG.err < /dev/null &
  PID=$!
  echo "{\"event\":\"launch\",\"attempt\":$attempt,\"pid\":$PID,\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  wait $PID; rc=$?
  if [ -f "$OUT/DONE_$STOP.json" ]; then echo "{\"event\":\"done\",\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 0; fi
  last=$($PY -c "import torch;print(torch.load('$OUT/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo -1)
  err=$(tail -n 3 $LOG.err | tr '\n"' '  ' | cut -c1-400)
  seg=$(dmesg 2>/dev/null | grep -c segfault)
  echo "{\"event\":\"interruption\",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"dmesg_segfaults_total\":$seg,\"stderr_tail\":\"$err\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  if [ "$last" = "$laststep" ]; then noprog=$((noprog+1)); else noprog=0; fi
  laststep=$last
  if [ $noprog -ge 3 ] && ! echo "$err" | grep -qiE 'out of memory|CUDA'; then echo "{\"event\":\"FATAL\",\"reason\":\"no progress in 4 attempts\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 9; fi
  pkill -P $$ -f a1_dp 2>/dev/null; sleep 30
done
echo "{\"event\":\"FATAL\",\"reason\":\"attempt budget exhausted\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 9
