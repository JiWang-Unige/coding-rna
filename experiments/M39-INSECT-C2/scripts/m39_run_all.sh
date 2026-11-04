#!/bin/bash
# M39 GPU sequence (protocol 29 §4.5, user-confirmed order): smoke (once) -> I10_k50U -> I10_k50F -> I10_T.  Each step is resumable and skipped when already done
# (SMOKE_OK / ARM_DONE); the sequence stops at the first failure (event in run_all.jsonl).  It does NOT wait for or check other GPU jobs: start it only from the
# coordinating session's queue (after k25U), e.g.
#   nohup bash ~/coding-rna-reorg/experiments/M39-INSECT-C2/scripts/m39_run_all.sh > ~/coding-rna/logs/m39-run-all.out 2>&1 < /dev/null &
# Status:  tail -n 5 ~/coding-rna/outputs/M39-INSECT-C2/run_all.jsonl; tail -n 3 ~/coding-rna/outputs/M39-INSECT-C2/I10_k50U/metrics.jsonl
cd $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts
source ./m39_env.sh
mkdir -p $OUTROOT
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $OUTROOT/run_all.jsonl; }
ev run_all_start
if [ ! -f $OUTROOT/smoke/SMOKE_OK ]; then
  bash m39_smoke.sh || { ev smoke_failed; exit 1; }
  ev smoke_ok
fi
for ARM in I10_k50U I10_k50F I10_T; do
  [ -f $OUTROOT/$ARM/ARM_DONE ] && continue
  ev arm_start ",\"arm\":\"$ARM\""
  bash m39_stage.sh $ARM; rc=$?
  [ -f $OUTROOT/$ARM/ARM_DONE ] || { ev arm_failed ",\"arm\":\"$ARM\",\"rc\":$rc"; exit 9; }
  ev arm_done ",\"arm\":\"$ARM\""
done
touch $OUTROOT/M39_GPU_DONE
ev run_all_done
