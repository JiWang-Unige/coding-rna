#!/bin/bash
# M40 GPU queue entry (waits for the M40_AFTER / M40_AFTER_ANY sentinels, then for a free GPU without TE):  GPU smoke -> M40_k25U_MLM -> M40_k50U_MLM -> M40_k0 -> Vipsania on the new-species pseudo-contigs (only without explicit arm list).
# Each step is resumable and retried up to 6 times after the GPU has >= 6 GB free three checks in a row (as m34_after_t2s3.sh); a smoke failure stops the queue.
# Prerequisites (CPU, done once): m40_repeats.py build ; m40_schedule.py.   Events: ~/coding-rna/outputs/M40-SEMISUP-MLM/queue.jsonl
# Options (env):
#   M40_AFTER="<f> <f> ..."      wait until ALL of these exist before anything starts (e.g. .../M39-INSECT-C2/M39_GPU_DONE)
#   M40_AFTER_ANY="<f> <f> ..."  then wait until ANY one of these exists (release on success OR on a failure marker of the previous queue).
#                                An entry 'file::regex' counts as present when `file` exists and contains a line matching `regex` (grep -E).
#   M34_MIN_FREE_MIB (default 6000)  free GPU memory needed for 3 checks 2 min apart; additionally no TE claim / TE process on the GPU (te_on_gpu, as M39/M41)
# usage: nohup bash m40_run_all.sh [arm ...] > /dev/null 2>&1 &
# Launch after M41 (user order 11-11: M39 -> M41 -> M40), on ont:
#   O41=$HOME/coding-rna/outputs/M41-INTRON-LENGTH; (M40_AFTER_ANY="$O41/M41_GPU_DONE $O41/SMOKE_FAIL $O41/smoke/SMOKE_FAIL $O41/queue.jsonl::gave_up" \
#     setsid nohup bash $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_run_all.sh > $HOME/coding-rna/logs/m40-run-all.out 2>&1 < /dev/null &)
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
mkdir -p $OUT40
LOG=$OUT40/queue.jsonl
exec 7>$OUT40/.queue.lock
flock -n 7 || { echo "M40 queue already running"; exit 3; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ARMS=${@:-M40_k25U_MLM M40_k50U_MLM M40_k0}
te_on_gpu() {   # 0 if a TE-project process holds the GPU or TE claims it (same test as m39_after_k25u.sh / m41_run_all.sh)
  local p
  [ -e $HOME/ab-initio-TE/outputs/TE_GPU_CLAIM ] && return 0
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE 'ab-initio-TE|te_benchmark|TE-GLM|/work/te' && return 0
    readlink /proc/$p/cwd 2>/dev/null | grep -qE 'ab-initio-TE|TE-GLM|/work/te' && return 0
    ps -o comm= -p $(ps -o ppid= -p $p 2>/dev/null) 2>/dev/null | grep -qiE 'apptainer|singularity' && return 0
  done
  return 1
}
gpu_wait() {
  local ok=0 n=0 free te
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if te_on_gpu; then te=1; else te=0; fi
    if [ $te -eq 0 ] && [ "${free:-0}" -ge ${M34_MIN_FREE_MIB:-6000} ]; then ok=$((ok+1)); else ok=0; fi
    [ $((n % 30)) -eq 0 ] && ev gpu_wait ",\"free_mib\":${free:-0},\"te\":$te,\"ok\":$ok"; n=$((n+1))
    [ $ok -lt 3 ] && sleep 120
  done
}
present() {   # file  or  file::regex
  case "$1" in
    *::*) [ -e "${1%%::*}" ] && grep -qE -- "${1#*::}" "${1%%::*}" ;;
    *) [ -e "$1" ] ;;
  esac
}
wait_after() {
  local f n hit
  for f in ${M40_AFTER:-}; do
    n=0; until present "$f"; do [ $((n % 30)) -eq 0 ] && ev waiting_after ",\"file\":\"$f\""; n=$((n+1)); sleep 120; done
  done
  if [ -n "${M40_AFTER_ANY:-}" ]; then
    n=0
    while :; do
      hit=; for f in $M40_AFTER_ANY; do present "$f" && { hit=$f; break; }; done
      [ -n "$hit" ] && break
      [ $((n % 30)) -eq 0 ] && ev waiting_after_any ",\"files\":\"$M40_AFTER_ANY\""; n=$((n+1)); sleep 120
    done
    ev after_any_released ",\"by\":\"$hit\""
  fi
  ev after_ok ",\"after\":\"${M40_AFTER:-}\",\"after_any\":\"${M40_AFTER_ANY:-}\""
}
step() {   # name sentinel command...
  local NAME=$1 SENT=$2 try rc; shift 2
  for try in 1 2 3 4 5 6; do
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$((try-1))"; return 0; }
    ev step_start ",\"step\":\"$NAME\",\"try\":$try,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
    "$@"; rc=$?
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$try"; return 0; }
    ev step_failed ",\"step\":\"$NAME\",\"try\":$try,\"rc\":$rc"
    [ $rc -eq 11 ] && break
    [ "$NAME" = smoke ] && [ -e $OUT40/smoke/SMOKE_FAIL ] && break     # a failed smoke check is not retried
    gpu_wait
  done
  ev step_gave_up ",\"step\":\"$NAME\""
  exit 9
}
ev queued ",\"order\":\"smoke $ARMS\""
wait_after
gpu_wait
step smoke $OUT40/smoke/SMOKE_OK bash $M40S/m40_smoke.sh
for A in $ARMS; do
  step $A $OUTROOT/pc2_$A/ARM_DONE bash $M40S/m40_stage.sh $A
done
# external 0 % point (decision m40-vipsania-newsp, user 2026-11-11): Vipsania pretrained-only on the 8 new-species pseudo-contigs; never blocks the queue
if [ -z "$1" ] && [ ! -e $OUT40/vip_newsp/SKIP ]; then
  ev step_start ",\"step\":\"vip_newsp\""
  bash $M40S/m40_vip_newsp.sh; rc=$?
  if [ -e $OUT40/vip_newsp/DONE ]; then ev step_done ",\"step\":\"vip_newsp\""; else ev step_failed ",\"step\":\"vip_newsp\",\"rc\":$rc"; fi
fi
ev all_done
