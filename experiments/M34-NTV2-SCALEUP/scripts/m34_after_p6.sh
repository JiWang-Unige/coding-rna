#!/bin/bash
# ont GPU queue after the PANEL-P6 baselines (user decision ont-order-after-k25u, 2026-11-09: P6 -> P4 seeds -> M39; decision p4-seed-robustness: one more seed
# for k50U and k50F, before M39).  M39 was stopped by hand at 14:38 on 11-09 (I10_k50U checkpoint step 320); it resumes here from that checkpoint.
#   0 wait      the ont P6 queue (p6_queue.sh) is gone AND PANEL-P6/P6_GPU_DONE_<host> exists (if the queue stops without it, keep waiting: the P6 session
#               restarts it; $O/AFTER_P6_GO_NOW skips this wait), then no TE process / claim on the GPU and >= 6 GB free for 3 checks
#   1 T10_k50U_s2   m34_stage.sh (seed 2: sampling_seed 2, lora_seed 2; same labeled genes as k50U)
#   2 T10_k50F_s2   m34_stage.sh (same)
#   3 M39       m39_run_all.sh (resumable: smoke done, I10_k50U from ckpt_latest.pt, then I10_k50F, I10_T)
# Events: $O/after_k50f.jsonl (same file as before, so the existing monitors keep working)
# Launch: cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && (setsid nohup bash m34_after_p6.sh > ~/coding-rna/logs/m34-after-p6.out 2>&1 < /dev/null &)
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
P6=$HOME/coding-rna/outputs/PANEL-P6
LOG=$O/after_k50f.jsonl
LOGS=$HOME/coding-rna/logs
exec 7>$O/.after_p6.lock
flock -n 7 || { echo "after-p6 queue already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
export M34_BACKEND=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}') M34_LORA=on
[ "$M34_BACKEND" = pc2 ] || { ev refused ",\"reason\":\"GO_STAGE1 backend is $M34_BACKEND\""; exit 1; }
source ./m34_env.sh
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1

te_on_gpu() {   # 0 if a TE-project process holds the GPU (same test as m39_after_k25u.sh)
  local p
  [ -e $HOME/ab-initio-TE/outputs/TE_GPU_CLAIM ] && return 0
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE 'ab-initio-TE|te_benchmark|TE-GLM|/work/te' && return 0
    readlink /proc/$p/cwd 2>/dev/null | grep -qE 'ab-initio-TE|TE-GLM|/work/te' && return 0
    ps -o comm= -p $(ps -o ppid= -p $p 2>/dev/null) 2>/dev/null | grep -qiE 'apptainer|singularity' && return 0
  done
  return 1
}

gpu_wait() {    # no TE on the GPU and >= M34_MIN_FREE_MIB free for 3 consecutive checks (10 min apart)
  local ok=0 n=0 free
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if ! te_on_gpu && [ "${free:-0}" -ge ${M34_MIN_FREE_MIB:-6000} ]; then ok=$((ok+1)); else ok=0; fi
    [ $((n % 6)) -eq 0 ] && [ $ok -eq 0 ] && ev waiting_gpu_free_of_te; n=$((n+1))
    sleep 600
  done
}

step() {   # name sentinel command...
  local NAME=$1 SENT=$2 try rc; shift 2
  for try in 1 2 3 4 5 6; do
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$((try-1))"; return 0; }
    gpu_wait
    ev step_start ",\"step\":\"$NAME\",\"try\":$try,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
    "$@"; rc=$?
    [ -e $SENT ] && { ev step_done ",\"step\":\"$NAME\",\"try\":$try"; return 0; }
    ev step_failed ",\"step\":\"$NAME\",\"try\":$try,\"rc\":$rc"
    [ $rc -eq 11 ] && break          # DEV refused / model changed: never retried
  done
  ev step_gave_up ",\"step\":\"$NAME\""
  exit 9
}

ev queued ",\"order\":\"wait P6 -> T10_k50U_s2 T10_k50F_s2 -> M39 (resume)\",\"decision\":\"ont-order-after-k25u + p4-seed-robustness\""
# 0 ---- wait for the ont P6 queue
n=0
until [ -e $O/AFTER_P6_GO_NOW ] || { ! pgrep -u $(id -un) -f 'p6_queue[.]sh' > /dev/null && [ -e $P6/P6_GPU_DONE_$(hostname) ]; }; do
  [ $((n % 30)) -eq 0 ] && ev waiting_p6 ",\"p6_queue_alive\":$(pgrep -u $(id -un) -f 'p6_queue[.]sh' > /dev/null && echo 1 || echo 0)"; n=$((n+1)); sleep 120
done
ev p6_wait_over ",\"p6_done\":$([ -e $P6/P6_GPU_DONE_$(hostname) ] && echo 1 || echo 0)"
# 1-2 ---- P4 seed replicates
step T10_k50U_s2 $O/pc2_T10_k50U_s2/ARM_DONE bash m34_stage.sh T10_k50U_s2
step T10_k50F_s2 $O/pc2_T10_k50F_s2/ARM_DONE bash m34_stage.sh T10_k50F_s2
# 3 ---- M39 (resumes I10_k50U from its checkpoint)
M=$HOME/coding-rna/outputs/M39-INSECT-C2
m39() { bash $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts/m39_run_all.sh >> $LOGS/m39-run-all.out 2>&1 < /dev/null; }
step M39 $M/M39_GPU_DONE m39
ev all_done
