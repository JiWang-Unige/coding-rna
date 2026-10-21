#!/bin/bash
# stage driver for one tau: train to 256/512/1024 (schedule always defined for 4608 steps), validate each model on the held-out
# TRAIN chromosomes (never DEV); early-stop the tau if the FIRST checkpoint (256) is below the gate.  usage: a2_stage.sh <tau> [stops...]
TAU=$1; shift; STOPS=${@:-"256 512 1024"}
O=$HOME/coding-rna/outputs/M29-A2-TEMPERATURE; R=$O/run_tau$TAU
cd $HOME/coding-rna-reorg/experiments/M29-A2-TEMPERATURE/scripts
for STOP in $STOPS; do
  bash a2_supervise.sh a2-train-tau$TAU $R $TAU 4608 $STOP 12
  V=$O/val_tau${TAU}_$STOP
  if [ ! -f $V/valscore.json ]; then
    bash a2_infer_retry.sh $V $R/model_$(printf %06d $STOP).pt NC_003071.7,NC_089044.1
    taskset -c 18-23 $HOME/envs/generanno/bin/python a2_valscore.py $V >> $O/valscore_tau$TAU.log 2>&1
  fi
  D=$(python3 -c "import json;print(json.load(open('$V/valscore.json'))['delta_F1'])")
  echo "{\"tau\":$TAU,\"stop\":$STOP,\"delta_F1_vs_gate\":$D,\"time\":\"$(date -Is)\"}" >> $O/stage_log.jsonl
  if [ "$STOP" = "256" ] && python3 -c "import sys;sys.exit(0 if $D<0 else 1)"; then
    echo "{\"tau\":$TAU,\"event\":\"early_stop_first_checkpoint_below_gate\"}" >> $O/stage_log.jsonl; break
  fi
done
