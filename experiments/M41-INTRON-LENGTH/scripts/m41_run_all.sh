#!/bin/bash
# M41 GPU entry (user decision m41-run = A, 2026-11-11: one seed of T2_ILS, r = 1, theta learnable; protocol 33).  Launched by the maintenance session
# (owner of the ont GPU queue) at the place it chooses in the queue; this script does not wait for other queue items itself, except the optional
# M41_AFTER sentinel list below, and it waits for a free GPU before each step.
#   1 smoke   m41_smoke.sh   -> outputs/M41-INTRON-LENGTH/smoke/SMOKE_OK   (failure stops here: no arm run on a failed smoke)
#   2 arm     m41_stage.sh   -> outputs/M41-INTRON-LENGTH/pc2_T2_ILS/ARM_DONE (resumable; DEV once; set B)
# Resumable: steps whose sentinel exists are skipped.  Events: ~/coding-rna/outputs/M41-INTRON-LENGTH/queue.jsonl
# Options (env):  M41_AFTER="<file> <file> ..."  wait until all these sentinels exist (e.g. .../M39_GPU_DONE or M40's done sentinel)
#                 M41_MIN_FREE_MB (default 30000) free GPU memory required for 3 checks 2 min apart, and no TE process / claim on the GPU
# Launch (ont):   cd ~/coding-rna-reorg/experiments/M41-INTRON-LENGTH/scripts && \
#                 (M41_AFTER="$HOME/coding-rna/outputs/M39-INSECT-C2/M39_GPU_DONE" setsid nohup bash m41_run_all.sh > ~/coding-rna/logs/m41-run-all.out 2>&1 < /dev/null &)
S41=$(cd "$(dirname "$0")" && pwd)
O=$HOME/coding-rna/outputs/M41-INTRON-LENGTH
LOG=$O/queue.jsonl
mkdir -p $O $HOME/coding-rna/logs
exec 7>$O/.run_all.lock
flock -n 7 || { echo "m41_run_all already running"; exit 3; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }

te_on_gpu() {
  local p
  [ -e $HOME/ab-initio-TE/outputs/TE_GPU_CLAIM ] && return 0
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE 'ab-initio-TE|te_benchmark|TE-GLM|/work/te' && return 0
    readlink /proc/$p/cwd 2>/dev/null | grep -qE 'ab-initio-TE|TE-GLM|/work/te' && return 0
  done
  return 1
}
gpu_wait() {
  local ok=0 n=0 free
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if ! te_on_gpu && [ "${free:-0}" -ge ${M41_MIN_FREE_MB:-30000} ]; then ok=$((ok+1)); else ok=0; fi
    [ $((n % 30)) -eq 0 ] && ev gpu_wait ",\"free_mb\":${free:-0},\"ok\":$ok"; n=$((n+1))
    [ $ok -lt 3 ] && sleep 120
  done
}

ev started ",\"decision\":\"m41-run=A\",\"after\":\"${M41_AFTER:-}\""
for f in ${M41_AFTER:-}; do
  n=0
  until [ -e "$f" ]; do [ $((n % 30)) -eq 0 ] && ev waiting_after ",\"file\":\"$f\""; n=$((n+1)); sleep 120; done
done
ev after_ok

if [ ! -e $O/smoke/SMOKE_OK ]; then
  gpu_wait; ev smoke_start
  bash $S41/m41_smoke.sh >> $O/smoke.out 2>&1 < /dev/null
  [ -e $O/smoke/SMOKE_OK ] || { ev smoke_failed; touch $O/SMOKE_FAIL; exit 1; }
fi
ev smoke_ok
if [ ! -e $O/pc2_T2_ILS/ARM_DONE ]; then
  for try in 1 2 3; do
    gpu_wait; ev arm_start ",\"try\":$try"
    bash $S41/m41_stage.sh >> $O/stage.out 2>&1 < /dev/null; rc=$?
    [ -e $O/pc2_T2_ILS/ARM_DONE ] && break
    ev arm_failed ",\"try\":$try,\"rc\":$rc"
    [ $rc -eq 11 ] && break
  done
fi
[ -e $O/pc2_T2_ILS/ARM_DONE ] && { touch $O/M41_GPU_DONE; ev all_done; } || { ev gave_up; exit 9; }
