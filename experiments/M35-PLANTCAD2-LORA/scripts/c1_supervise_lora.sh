#!/bin/bash
# usage: c1_supervise_lora.sh <name> <outdir> <total_steps> <stop_at> [extra c1_train_lora args...]
# crash-resilient trainer launcher (M35): resumes from ckpt_latest.pt after any interruption (GPU OOM from other users, segfault, ...) and logs each one to interruptions.jsonl.
NAME=$1; OUT=$2; STEPS=$3; STOP=$4; shift 4
CPUS=${M35_CPUS:-2,4,6,8,10,12,14,16-24}
WORKERS=${M35_WORKERS:-12}
MINFREE=${M35_MIN_FREE_MIB:-9000}
INIT=$HOME/coding-rna/outputs/M32-C1-PLANTCAD2/runA/model_002048.pt
PY=$HOME/envs/plantcad2/bin/python
LOG=$HOME/coding-rna/logs/$NAME-$(date +%Y%m%d-%H%M%S)
mkdir -p "$OUT"
cd $HOME/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA/scripts
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
for attempt in $(seq 1 60); do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt "$MINFREE" ]; then echo "{\"event\":\"gpu_busy_wait\",\"free_mib\":$free,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; sleep 120; continue; fi
  taskset -c $CPUS $PY c1_train_lora.py run --out "$OUT" --init $INIT --steps $STEPS --stop-at $STOP --batch 8 --workers $WORKERS --cpus $CPUS --resume "$@" >> $LOG.out 2>> $LOG.err < /dev/null &
  PID=$!
  echo "{\"event\":\"launch\",\"attempt\":$attempt,\"pid\":$PID,\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  wait $PID; rc=$?
  if [ -f "$OUT/DONE_$STOP.json" ]; then echo "{\"event\":\"done\",\"stop_at\":$STOP,\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"; break; fi
  last=$($PY -c "import torch;print(torch.load('$OUT/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo none)
  err=$(tail -n 3 $LOG.err | tr '\n"' '  ' | cut -c1-400)
  seg=$(dmesg 2>/dev/null | grep -c segfault)
  echo "{\"event\":\"interruption\",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"dmesg_segfaults_total\":$seg,\"stderr_tail\":\"$err\",\"time\":\"$(date -Is)\"}" >> "$OUT/interruptions.jsonl"
  pkill -P $$ -f a1_dp 2>/dev/null; sleep 30
done
