#!/bin/bash
# stage driver: train to each stop (schedule always defined for STEPS), validate on the held-out TRAIN chromosomes (never DEV), early-stop on patience.
# usage: c1_stage.sh <runname> <steps> <patience> <layer> [first_stop] [stop_every]
RUN=$1; STEPS=${2:-4608}; PATIENCE=${3:-2}; LAYER=${4:--1}; FIRST=${5:-512}; EVERY=${6:-512}
O=$HOME/coding-rna/outputs/M32-C1-PLANTCAD2; R=$O/$RUN
REGIONS=NC_003071.7,NC_089044.1        # full held-out TRAIN chromosomes (Arabidopsis chr2, rice NC_089044.1): internal validation, never DEV
PY=$HOME/envs/plantcad2/bin/python
cd $HOME/coding-rna-reorg/experiments/M32-C1-PLANTCAD2/scripts
best=-1; bad=0
for STOP in $(seq $FIRST $EVERY $STEPS); do
  bash c1_supervise.sh m32-train-$RUN $R $STEPS $STOP --layer $LAYER
  V=$R/val_$STOP
  if [ ! -f $V/valscore.json ]; then
    bash c1_infer_retry.sh $V $R/model_$(printf %06d $STOP).pt $REGIONS 4
    taskset -c ${M32_SCORE_CPUS:-18-23} $PY c1_val.py $V $STOP $REGIONS --metrics $R/metrics.jsonl >> $R/val.log 2>&1
  fi
  F=$(python3 -c "import json;print(json.load(open('$V/valscore.json')).get('pooled_F1',0))")
  echo "{\"stop\":$STOP,\"val_pooled_F1\":$F,\"time\":\"$(date -Is)\"}" >> $R/stage_log.jsonl
  if python3 -c "import sys;sys.exit(0 if $F>$best else 1)"; then best=$F; bad=0; echo $STOP > $R/best_step.txt; else bad=$((bad+1)); fi
  rm -rf $V/pot
  if [ $bad -ge $PATIENCE ]; then echo "{\"event\":\"early_stop\",\"stop\":$STOP,\"best\":$best,\"time\":\"$(date -Is)\"}" >> $R/stage_log.jsonl; break; fi
done
echo "{\"event\":\"stage_complete\",\"best_step\":\"$(cat $R/best_step.txt)\",\"time\":\"$(date -Is)\"}" >> $R/stage_log.jsonl
