#!/bin/bash
# ont GPU queue after the PANEL-P6 baselines, revised 2026-11-17 (user: "均按推荐来" on decisions c1-fail-next and insect-lp-slot):
#   P4 seed replicates -> M41 (intron-length control, moved forward after C1 failed on the multi-exon deficit) -> M39 NT-v2 three arms (resume)
#   -> protocol 32 arm L (one-hot floor) -> arm P (PlantCAD2) -> M40 (semi-supervised MLM).
# Replaces m34_after_p6.sh + the separately launched m41_run_all.sh / m40_run_all.sh waiters (stopped by hand while all were only waiting).
#   1-2 T10_k50U_s2, T10_k50F_s2   m34_stage.sh (pc2), retried like m34_after_p6.sh
#   3   M41                         m41_run_all.sh (own GPU / TE gate; a smoke failure stops M41 but not this queue)
#   4   M39                         m39_run_all.sh (NT-v2; resumes I10_k50U from its checkpoint)
#   5-6 arm L, arm P                M39_BACKEND=onehot / pc2  m39_stage.sh I10_k50U
#   7   M40                         m40_run_all.sh (own GPU / TE gate)
# Every step waits for: no TE claim / TE process on the GPU and >= M34_MIN_FREE_MIB free for 3 checks (10 min apart).
# Events: $O/after_k50f.jsonl (same file as before, so the existing monitors keep working)
# Launch: cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && (setsid nohup bash m34_after_p6b.sh > ~/coding-rna/logs/m34-after-p6b.out 2>&1 < /dev/null &)
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
LOG=$O/after_k50f.jsonl
LOGS=$HOME/coding-rna/logs
EXP=$HOME/coding-rna-reorg/experiments
exec 7>$O/.after_p6.lock
flock -n 7 || { echo "after-p6 queue already running"; exit 3; }
cd $EXP/M34-NTV2-SCALEUP/scripts
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
export M34_BACKEND=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}') M34_LORA=on
[ "$M34_BACKEND" = pc2 ] || { ev refused ",\"reason\":\"GO_STAGE1 backend is $M34_BACKEND\""; exit 1; }
source ./m34_env.sh
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1

te_on_gpu() {   # 0 if a TE-project process holds the GPU (same test as m39_after_k25u.sh); an empty claim file is not a claim (as the P6 gate, 134926f)
  local p
  [ -s $HOME/ab-initio-TE/outputs/TE_GPU_CLAIM ] && return 0
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE 'ab-initio-TE|te_benchmark|TE-GLM|/work/te' && return 0
    readlink /proc/$p/cwd 2>/dev/null | grep -qE 'ab-initio-TE|TE-GLM|/work/te' && return 0
    ps -o comm= -p $(ps -o ppid= -p $p 2>/dev/null) 2>/dev/null | grep -qiE 'apptainer|singularity' && return 0
  done
  return 1
}

gpu_wait() {
  local ok=0 n=0 free
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if ! te_on_gpu && [ "${free:-0}" -ge ${M34_MIN_FREE_MIB:-6000} ]; then ok=$((ok+1)); else ok=0; fi
    [ $((n % 6)) -eq 0 ] && [ $ok -eq 0 ] && ev waiting_gpu_free_of_te; n=$((n+1))
    sleep 600
  done
}

step() {   # name sentinel command...   (blocking: gives up the whole queue after 6 tries)
  local NAME=$1 SENT=$2 try rc; shift 2
  for try in 1 2 3 4 5 6; do
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$((try-1))"; return 0; }
    gpu_wait
    ev step_start ",\"step\":\"$NAME\",\"try\":$try,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
    "$@"; rc=$?
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$try"; return 0; }
    ev step_failed ",\"step\":\"$NAME\",\"try\":$try,\"rc\":$rc"
    [ $rc -eq 11 ] && break
  done
  ev step_gave_up ",\"step\":\"$NAME\""
  exit 9
}

soft() {   # name command...   (non-blocking: the sub-queue has its own gate and retries; its failure does not stop this queue)
  local NAME=$1 rc; shift
  ev step_start ",\"step\":\"$NAME\""
  "$@"; rc=$?
  ev step_end ",\"step\":\"$NAME\",\"rc\":$rc"
}

m39arm() { ( cd $EXP/M39-INSECT-C2/scripts && M39_BACKEND=$1 bash m39_stage.sh I10_k50U >> $LOGS/m39-$1-I10_k50U.out 2>&1 < /dev/null ); }
m39all() { ( cd $EXP/M39-INSECT-C2/scripts && M39_BACKEND=ntv2 bash m39_run_all.sh >> $LOGS/m39-run-all.out 2>&1 < /dev/null ); }
m41all() { ( cd $EXP/M41-INTRON-LENGTH/scripts && M41_MIN_FREE_MB=6000 bash m41_run_all.sh >> $LOGS/m41-run-all.out 2>&1 < /dev/null ); }
m40all() { bash $EXP/M40-SEMISUP-MLM/scripts/m40_run_all.sh >> $LOGS/m40-run-all.out 2>&1 < /dev/null; }

M39O=$HOME/coding-rna/outputs/M39-INSECT-C2
ev queued ",\"order\":\"T10_k50U_s2 T10_k50F_s2 -> M41 -> M39 (resume) -> L -> P -> M40\",\"decision\":\"c1-fail-next + insect-lp-slot (11-17)\",\"replaces\":\"m34_after_p6.sh, m41_run_all.sh waiter, m40_run_all.sh waiter\""
step T10_k50U_s2 $O/pc2_T10_k50U_s2/ARM_DONE bash m34_stage.sh T10_k50U_s2
step T10_k50F_s2 $O/pc2_T10_k50F_s2/ARM_DONE bash m34_stage.sh T10_k50F_s2
gpu_wait; soft M41 m41all
step M39 $M39O/M39_GPU_DONE m39all
step M39_L $M39O/oh_I10_k50U/ARM_DONE m39arm onehot
step M39_P $M39O/pc2_I10_k50U/ARM_DONE m39arm pc2
gpu_wait; soft M40 m40all
ev all_done
