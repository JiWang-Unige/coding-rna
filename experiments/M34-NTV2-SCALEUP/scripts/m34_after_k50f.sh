#!/bin/bash
# ont GPU queue after M34 k50F (user decision gpu-after-k50f, 2026-11-01: k50F -> panel inference with the final model T2 -> +2 seeds of T2 -> k25U).  It replaces the
# k25U half of m34_after_m36.sh: when this queue is started, that chain process and m34_relaunch_after_m36.sh are stopped by hand (the running k50F stage keeps running).
# Steps run one after another, each resumable through its own sentinel; a failing step (e.g. m34_supervise.sh giving up while other users hold the GPU) is retried after
# the GPU has >= 7.5 GB free for 3 consecutive checks (120 s apart), at most 6 tries per step; then the queue stops.
#   1 k50F   wait until the k50F stage process is gone; m34_stage.sh T10_k50F again if pc2_T10_k50F/ARM_DONE is missing (resumes)
#   2 panel  m34_panel.sh (P5, core 10 species, T2 step 2304)                                      -> panel/PANEL_INFER_DONE
#   3 seeds  m34_stage.sh T2_s2, zero-shot set-B inference of its selected checkpoint, then T2_s3 the same   -> pc2_T2_sN/ARM_DONE, pc2_T2_sN/setB/SETB_DONE
#   4 k25U   m34_stage.sh T10_k25U                                                                 -> pc2_T10_k25U/ARM_DONE
# Events: $O/after_k50f.jsonl
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
LOG=$O/after_k50f.jsonl
LOGS=$HOME/coding-rna/logs
exec 7>$O/.after_k50f.lock
flock -n 7 || { echo "after-k50f queue already running"; exit 3; }
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

ev queued ",\"order\":\"T10_k50F panel T2_s2 T2_s3 T10_k25U\""
# 1 ---- k50F: let the running stage finish, re-run it (resumable) only if it stopped before ARM_DONE
n=0
while pgrep -u $(id -un) -f 'm34_stage.sh T10_k50[F]' > /dev/null; do
  [ $((n % 30)) -eq 0 ] && ev waiting_k50f; n=$((n+1)); sleep 120
done
step T10_k50F $O/pc2_T10_k50F/ARM_DONE bash m34_stage.sh T10_k50F
# 2 ---- panel inference (P5)
step panel $O/panel/PANEL_INFER_DONE bash m34_panel.sh
# 3 ---- +2 seeds of T2
for ARM in T2_s2 T2_s3; do
  [ -f $O/schedules/schedule_$ARM.jsonl ] || { ev refused ",\"reason\":\"schedule_$ARM.jsonl missing\""; exit 1; }
  step $ARM $O/pc2_$ARM/ARM_DONE bash m34_stage.sh $ARM
  step ${ARM}_setB $O/pc2_$ARM/setB/SETB_DONE setb $ARM
done
# 4 ---- k25U
step T10_k25U $O/pc2_T10_k25U/ARM_DONE bash m34_stage.sh T10_k25U
ev all_done
