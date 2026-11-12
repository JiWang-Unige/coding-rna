#!/bin/bash
# usage: m40_supervise.sh <arm> <outdir> <stop_at> [extra m40_train args...]
# = m34_supervise.sh with the M40 trainer: resumes from ckpt_latest.pt after any interruption, waits while the GPU has < 6.0 GB free, FATAL after 4 attempts without progress.
ARM=$1; OUT=$2; STOP=$3; shift 3
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
CPUS=${M34_CPUS:-2,4,6,10,12,14,16-24}
LOG=$HOME/coding-rna/logs/m40-train-$ARM-$(date +%Y%m%d-%H%M%S)
mkdir -p "$OUT" $HOME/coding-rna/logs
cd $M40S
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
laststep=-1; noprog=0
for attempt in $(seq 1 80); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt ${M34_MIN_FREE_MIB:-6000} ]; then echo "{\"event\":\"gpu_busy_wait\",\"free_mib\":$free,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; sleep 120; continue; fi
  taskset -c $CPUS $PY m40_train.py run --arm $ARM --out "$OUT" --steps 3072 --stop-at $STOP --batch 8 --workers 12 --cpus $CPUS --resume "$@" >> $LOG.out 2>> $LOG.err < /dev/null &
  PID=$!
  echo "{\"event\":\"launch\",\"attempt\":$attempt,\"pid\":$PID,\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  wait $PID; rc=$?
  if [ -f "$OUT/DONE_$STOP.json" ]; then echo "{\"event\":\"done\",\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 0; fi
  last=$($PY -c "import torch;print(torch.load('$OUT/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo -1)
  err=$(tail -n 3 $LOG.err | tr '\n"' '  ' | cut -c1-400)
  echo "{\"event\":\"interruption\",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"stderr_tail\":\"$err\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  if [ "$last" = "$laststep" ]; then noprog=$((noprog+1)); else noprog=0; fi
  laststep=$last
  if [ $noprog -ge 3 ] && ! echo "$err" | grep -qiE 'out of memory|CUDA'; then echo "{\"event\":\"FATAL\",\"reason\":\"no progress in 4 attempts\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 9; fi
  pkill -P $$ -f a1_dp 2>/dev/null; sleep 30
done
echo "{\"event\":\"FATAL\",\"reason\":\"attempt budget exhausted\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; exit 9
