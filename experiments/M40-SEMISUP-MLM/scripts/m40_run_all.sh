#!/bin/bash
# M40 GPU queue entry (the maintainer starts it after M39; it does NOT wait for M39 itself):  GPU smoke -> M40_k25U_MLM -> M40_k50U_MLM -> M40_k0 -> Vipsania on the new-species pseudo-contigs (only without explicit arm list).
# Each step is resumable and retried up to 6 times after the GPU has >= 6 GB free three checks in a row (as m34_after_t2s3.sh); a smoke failure stops the queue.
# Prerequisites (CPU, done once): m40_repeats.py build ; m40_schedule.py.   Events: ~/coding-rna/outputs/M40-SEMISUP-MLM/queue.jsonl
# usage: nohup bash m40_run_all.sh [arm ...] > /dev/null 2>&1 &
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
mkdir -p $OUT40
LOG=$OUT40/queue.jsonl
exec 7>$OUT40/.queue.lock
flock -n 7 || { echo "M40 queue already running"; exit 3; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ARMS=${@:-M40_k25U_MLM M40_k50U_MLM M40_k0}
gpu_wait() {
  local ok=0 free
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "${free:-0}" -ge ${M34_MIN_FREE_MIB:-6000} ]; then ok=$((ok+1)); else ok=0; fi
    sleep 120
  done
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
