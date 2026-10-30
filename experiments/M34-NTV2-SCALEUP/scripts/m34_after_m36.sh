#!/bin/bash
# Run the remaining M34 arms T10_k50F then T10_k25U (pc2) after the M36 queue has finished (user decisions m34-remaining-arms 2026-10-27 and m34-k-workers 2026-10-30:
# 12 CRF workers at default priority, same as k50U).  Start condition (checked every 120 s): M36_DONE exists and no coding-rna GPU process for 3 consecutive checks,
# so the two coding-rna GPU jobs never overlap.  Arms run one after another; a failing arm stops the chain.
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
M36=$HOME/coding-rna/outputs/M36-BOOTSTRAP-ADAPT
LOG=$O/after_m36.jsonl
exec 7>$O/.after_m36.lock
flock -n 7 || { echo "after-m36 queue already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
source $HOME/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts/m36_gate.sh      # m36_gpu_blockers (coding-rna GPU pids of this user)
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ev queued ",\"arms\":\"T10_k50F T10_k25U\""
ok=0; n=0
while :; do
  if [ -f $M36/M36_DONE ] && [ -z "$(m36_gpu_blockers)" ]; then ok=$((ok+1)); else ok=0; fi
  [ $ok -ge 3 ] && break
  [ $((n % 30)) -eq 0 ] && ev waiting ",\"ok\":$ok"
  n=$((n+1)); sleep 120
done
export M34_BACKEND=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}') M34_LORA=on
source ./m34_env.sh
for ARM in T10_k50F T10_k25U; do
  ev arm_start ",\"arm\":\"$ARM\",\"backend\":\"$M34_BACKEND\",\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
  bash m34_stage.sh $ARM; rc=$?
  ev arm_end ",\"arm\":\"$ARM\",\"rc\":$rc"
  [ $rc -eq 0 ] || exit $rc
done
ev all_done
