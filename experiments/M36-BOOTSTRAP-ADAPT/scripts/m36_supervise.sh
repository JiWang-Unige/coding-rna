#!/bin/bash
# usage: m36_supervise.sh <arm> <species> <outdir> <stop_at> [extra m36_train.py args...]
# Crash-resilient launcher of one M36 adaptation run (the m34_supervise.sh pattern): resumes from ckpt_latest.pt after any interruption (GPU OOM caused by other users,
# segfault, ...), logs each one in <outdir>/interruptions.jsonl, waits while the GPU has < 7.5 GB free, and gives up (FATAL, exit 9) when 4 consecutive attempts make no
# progress for a reason other than CUDA / out-of-memory (a deterministic code error).  Exit 11 = m36_train.py refused (fingerprint of seeds / base / schedule changed).
ARM=$1; SP=$2; OUT=$3; STOP=$4; shift 4
export M34_BACKEND=pc2
PY=${M36_PY_TRAIN:-$HOME/envs/plantcad2/bin/python}
TRAIN=${M36_TRAIN_SCRIPT:-m36_train.py}
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
LOGS=${M36_LOGS:-$HOME/coding-rna/logs}
LOG=$LOGS/m36-train-${ARM}_${SP}-$(date +%Y%m%d-%H%M%S)
SLEEP_BUSY=${M36_SLEEP_BUSY:-120}; SLEEP_RETRY=${M36_SLEEP_RETRY:-30}
mkdir -p "$OUT" $LOGS
cd ${M36_SCRIPTS:-$HOME/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts}
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
free_mib() { if [ -n "$M36_FAKE_FREE_MIB" ]; then echo $M36_FAKE_FREE_MIB; else nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; fi; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> "$OUT/interruptions.jsonl"; }
laststep=-1; noprog=0
for attempt in $(seq 1 80); do
  free=$(free_mib)
  if [ "${free:-0}" -lt 7500 ]; then ev gpu_busy_wait ",\"free_mib\":$free"; sleep $SLEEP_BUSY; continue; fi
  taskset -c $CPUS $PY $TRAIN run --arm $ARM --species $SP --out "$OUT" --stop-at $STOP --workers 12 --cpus $CPUS --resume "$@" >> $LOG.out 2>> $LOG.err < /dev/null &
  PID=$!
  ev launch ",\"attempt\":$attempt,\"pid\":$PID,\"stop_at\":$STOP"
  wait $PID; rc=$?
  if [ -f "$OUT/DONE_$STOP.json" ]; then ev done ",\"stop_at\":$STOP"; exit 0; fi
  if [ $rc -eq 11 ]; then ev refused ",\"rc\":11"; exit 11; fi
  last=$($PY -c "import torch;print(torch.load('$OUT/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo -1)
  err=$(tail -n 3 $LOG.err | tr '\n"\\' '   ' | cut -c1-400)
  seg=$(dmesg 2>/dev/null | grep -c segfault)
  ev interruption ",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"dmesg_segfaults_total\":$seg,\"stderr_tail\":\"$err\""
  if [ "$last" = "$laststep" ]; then noprog=$((noprog+1)); else noprog=0; fi
  laststep=$last
  if [ $noprog -ge 3 ] && ! echo "$err" | grep -qiE 'out of memory|CUDA'; then ev FATAL ",\"reason\":\"no progress in 4 attempts\""; exit 9; fi
  pkill -P $$ -f a1_dp 2>/dev/null; sleep $SLEEP_RETRY
done
ev FATAL ",\"reason\":\"attempt budget exhausted\""; exit 9
