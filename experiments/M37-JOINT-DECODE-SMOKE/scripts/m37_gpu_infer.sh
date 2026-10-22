#!/bin/bash
# GPU stage of M37: model B inference on the 9311 chromosome 10 (CM000135.2) only (A = NC_089044.1 potentials are reused from M33 val).
# Protocol §2 constraints: start only if free GPU memory >= 8 GB; GPU phase <= 20 min (watchdog kills the process group once 1200 s pass
# without both potential files written); 4 decode workers on cores 25-31 (core 24 is in M35's affinity set), nice 10.  M34/M35 untouched.
O=$HOME/coding-rna/outputs/M37-JOINT-DECODE-SMOKE/inferB
MODEL=$HOME/coding-rna/outputs/M31-C1-NTV2/runB/model_001536.pt
CPUS=25,26,27,28,29,30,31
mkdir -p $O $HOME/coding-rna/logs
cd $HOME/coding-rna-reorg/experiments/M37-JOINT-DECODE-SMOKE/scripts
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
echo "{\"event\":\"gpu_check\",\"free_MiB\":$free,\"time\":\"$(date -Is)\"}" >> $O/gpu_log.jsonl
if [ "$free" -lt 8192 ]; then echo "{\"event\":\"abort_low_gpu_memory\",\"free_MiB\":$free}" >> $O/gpu_log.jsonl; exit 2; fi
setsid nice -n 10 taskset -c $CPUS $HOME/envs/generanno/bin/python m37_infer.py --model $MODEL --regions CM000135.2 --out $O --decode-workers 4 --cpus $CPUS > $HOME/coding-rna/logs/m37-inferB.out 2> $HOME/coding-rna/logs/m37-inferB.err < /dev/null &
PID=$!
START=$(date +%s); GPU_DONE=""
echo "{\"event\":\"started\",\"pid\":$PID,\"time\":\"$(date -Is)\"}" >> $O/gpu_log.jsonl
while kill -0 $PID 2>/dev/null; do
  sleep 20
  if [ -z "$GPU_DONE" ] && [ -f $O/pot/CM000135.2_+.npy ] && [ -f $O/pot/CM000135.2_-.npy ]; then GPU_DONE=1; echo "{\"event\":\"gpu_phase_done\",\"elapsed_s\":$(( $(date +%s)-START )),\"time\":\"$(date -Is)\"}" >> $O/gpu_log.jsonl; fi
  if [ -z "$GPU_DONE" ] && [ $(( $(date +%s)-START )) -gt 1200 ]; then
    echo "{\"event\":\"watchdog_kill_20min\",\"time\":\"$(date -Is)\"}" >> $O/gpu_log.jsonl
    kill -TERM -- -$PID; sleep 5; kill -KILL -- -$PID 2>/dev/null; exit 3
  fi
done
wait $PID; rc=$?
echo "{\"event\":\"exit\",\"rc\":$rc,\"elapsed_s\":$(( $(date +%s)-START )),\"time\":\"$(date -Is)\"}" >> $O/gpu_log.jsonl
