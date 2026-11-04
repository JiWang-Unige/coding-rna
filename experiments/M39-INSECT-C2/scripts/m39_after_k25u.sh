#!/bin/bash
# ont GPU queue continuation for M39-INSECT-C2 (user decisions animal-c2 + animal-c2-protocol, 2026-11-04: insects run after k25U).  The M34 queue
# m34_after_t2s2.sh (T2_s2 -> Vipsania -> T2_s3 -> k25U) is not modified; this script waits for it and then runs m39_run_all.sh
# (smoke -> I10_k50U -> I10_k50F -> I10_T).
#   start condition: pc2_T10_k25U/ARM_DONE exists and the M34 queue process is gone.  If the M34 queue ends without k25U done, nothing is started
#   (event m39_not_started) -- the coordinating session decides.
# Events: the M34 queue file after_k50f.jsonl (so the existing monitors see them) and M39's own run_all.jsonl.
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
LOG=$O/after_k50f.jsonl
exec 7>$HOME/coding-rna/outputs/.m39_after_k25u.lock
flock -n 7 || { echo "m39 continuation already running"; exit 3; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ev m39_queued ",\"after\":\"T10_k25U\""
n=0
while pgrep -u $(id -un) -f 'm34_after_t2s2[.]sh' > /dev/null; do
  [ $((n % 30)) -eq 0 ] && ev waiting_m34_queue; n=$((n+1)); sleep 120
done
[ -e $O/pc2_T10_k25U/ARM_DONE ] || { ev m39_not_started ",\"reason\":\"M34 queue ended without pc2_T10_k25U/ARM_DONE\""; exit 1; }
ok=0                                       # GPU free >= 7.5 GB for 3 consecutive checks (120 s apart), as in the M34 queue
while [ $ok -lt 3 ]; do
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "${free:-0}" -ge 7500 ]; then ok=$((ok+1)); else ok=0; fi
  sleep 120
done
ev step_start ",\"step\":\"M39\",\"try\":1,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
bash $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts/m39_run_all.sh >> $HOME/coding-rna/logs/m39-run-all.out 2>&1 < /dev/null; rc=$?
if [ -e $HOME/coding-rna/outputs/M39-INSECT-C2/M39_GPU_DONE ]; then ev step_done ",\"step\":\"M39\",\"try\":1"; else ev step_failed ",\"step\":\"M39\",\"try\":1,\"rc\":$rc"; fi
