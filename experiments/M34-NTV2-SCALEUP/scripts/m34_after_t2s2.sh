#!/bin/bash
# ont GPU queue from T2_s2 on (user decision vipsania-run, 2026-11-04: Vipsania on DEV + set B goes between T2_s2 and T2_s3).  It replaces m34_after_k50f.sh, whose
# process is stopped by hand while the T2_s2 stage keeps running; the remaining steps are the same, with the Vipsania step inserted:
#   1 T2_s2     wait until the T2_s2 stage process is gone; m34_stage.sh T2_s2 again if ARM_DONE is missing (resumes); then its zero-shot set-B inference
#   2 vipsania  wait for M38-VIPSANIA/READY (written by hand once the CPU preparation is verified), then vip_devb.sh -> M38-VIPSANIA/DEVB_DONE.
#               A failing or skipped Vipsania step (M38-VIPSANIA/SKIP) does not block the queue.
#   3 T2_s3     m34_stage.sh T2_s3, then its zero-shot set-B inference
#   4 k25U      m34_stage.sh T10_k25U
# Events: $O/after_k50f.jsonl (same file as before, so the existing monitors keep working)
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
V=$HOME/coding-rna/outputs/M38-VIPSANIA
VS=$HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts
LOG=$O/after_k50f.jsonl
LOGS=$HOME/coding-rna/logs
exec 7>$O/.after_t2s2.lock
flock -n 7 || { echo "after-t2s2 queue already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
export M34_BACKEND=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}') M34_LORA=on
[ "$M34_BACKEND" = pc2 ] || { ev refused ",\"reason\":\"GO_STAGE1 backend is $M34_BACKEND\""; exit 1; }
source ./m34_env.sh
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
REGB=NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1,NW_028238309.1,NW_028238461.1,NW_028238483.1     # evaluation set B (13_M30_protocol section 2), as P2b

gpu_wait() {
  local ok=0 free
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "${free:-0}" -ge 7500 ]; then ok=$((ok+1)); else ok=0; fi
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
    [ $rc -eq 11 ] && break          # DEV refused / model changed: never retried
    gpu_wait
  done
  ev step_gave_up ",\"step\":\"$NAME\""
  exit 9
}

setb() {   # arm : zero-shot set-B inference of the arm's selected checkpoint (same path as P2b; potentials not kept)
  local ARM=$1 D=$O/pc2_$1/setB B a free
  B=$(head -c 20 $O/pc2_$ARM/best_step.txt 2>/dev/null | tr -dc 0-9)
  [ -n "$B" ] || { ev setb_refused ",\"arm\":\"$ARM\",\"reason\":\"no best_step.txt\""; return 11; }
  mkdir -p $D
  for a in $(seq 1 60); do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "${free:-0}" -lt 5000 ]; then sleep 120; continue; fi
    taskset -c $CPUS $PY m34_p2b_infer.py --model $O/pc2_$ARM/model_$(printf %06d $B).pt --regions $REGB --out $D --decode-workers 4 --cpus $CPUS \
      >> $LOGS/m34-setB-$ARM.out 2>> $LOGS/m34-setB-$ARM.err < /dev/null
    if ls $D/DONE_*.json >/dev/null 2>&1; then echo "{\"arm\":\"$ARM\",\"step\":$B,\"time\":\"$(date -Is)\"}" > $D/SETB_DONE; return 0; fi
    sleep 30
  done
  return 9
}

ev queued ",\"order\":\"T2_s2 vipsania T2_s3 T10_k25U\",\"replaces\":\"m34_after_k50f.sh\""
# 1 ---- T2_s2: let the running stage finish, re-run it (resumable) only if it stopped before ARM_DONE; then its set-B inference
n=0
while pgrep -u $(id -un) -f 'm34_stage.sh T2_s[2]' > /dev/null; do
  [ $((n % 30)) -eq 0 ] && ev waiting_t2s2; n=$((n+1)); sleep 120
done
step T2_s2 $O/pc2_T2_s2/ARM_DONE bash m34_stage.sh T2_s2
step T2_s2_setB $O/pc2_T2_s2/setB/SETB_DONE setb T2_s2
# 2 ---- Vipsania on DEV + set B (M38); never blocks the queue
mkdir -p $V
n=0
while [ ! -e $V/READY ] && [ ! -e $V/SKIP ]; do
  [ $((n % 12)) -eq 0 ] && ev waiting_vipsania_ready; n=$((n+1)); sleep 300
done
if [ -e $V/SKIP ]; then
  ev step_skipped ",\"step\":\"vipsania\""
elif [ ! -e $V/DEVB_DONE ]; then
  ev step_start ",\"step\":\"vipsania\",\"try\":1,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
  bash $VS/vip_devb.sh >> $LOGS/m38-vipsania-devb.out 2>> $LOGS/m38-vipsania-devb.err < /dev/null; rc=$?
  if [ -e $V/DEVB_DONE ]; then ev step_done ",\"step\":\"vipsania\",\"try\":1"; else ev step_failed ",\"step\":\"vipsania\",\"try\":1,\"rc\":$rc"; fi
else
  ev step_done ",\"step\":\"vipsania\",\"try\":0"
fi
# 3 ---- T2_s3
ARM=T2_s3
[ -f $O/schedules/schedule_$ARM.jsonl ] || { ev refused ",\"reason\":\"schedule_$ARM.jsonl missing\""; exit 1; }
step $ARM $O/pc2_$ARM/ARM_DONE bash m34_stage.sh $ARM
step ${ARM}_setB $O/pc2_$ARM/setB/SETB_DONE setb $ARM
# 4 ---- k25U
step T10_k25U $O/pc2_T10_k25U/ARM_DONE bash m34_stage.sh T10_k25U
ev all_done
