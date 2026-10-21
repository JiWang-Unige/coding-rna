#!/bin/bash
# usage: m34_supervise.sh <arm> <outdir> <stop_at> [extra m34_train_lora args...]
# crash-resilient trainer launcher (as M31 c1_supervise_lora.sh): resumes from ckpt_latest.pt after any interruption (GPU OOM from other users, segfault, ...), logs each one in
# interruptions.jsonl, waits while the GPU has < 7.5 GB free, and gives up (FATAL) when 4 consecutive attempts make no progress (a deterministic code error).
ARM=$1; OUT=$2; STOP=$3; shift 3
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
LOG=$HOME/coding-rna/logs/m34-train-$ARMPFX$ARM-$(date +%Y%m%d-%H%M%S)
mkdir -p "$OUT" $HOME/coding-rna/logs
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
source ./m34_env.sh
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
laststep=-1; noprog=0
for attempt in $(seq 1 80); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt 7500 ]; then echo "{\"event\":\"gpu_busy_wait\",\"free_mib\":$free,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; sleep 120; continue; fi
  taskset -c $CPUS $PY m34_train_lora.py run --arm $ARM --out "$OUT" --steps 3072 --stop-at $STOP --batch 8 --workers 12 --cpus $CPUS --resume ${M34_LORA:+--lora $M34_LORA} "$@" >> $LOG.out 2>> $LOG.err < /dev/null &
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
