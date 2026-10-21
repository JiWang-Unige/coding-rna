#!/bin/bash
# Queue driver (runs detached on ont): wait until NO process whose command line contains M31-C1-NTV2 is left (M31 trainer / inference / supervisor / stage), then
# start the M32 stage (train -> internal validation -> early stop) exactly once.  One GPU job of ours at a time.
# usage: start_after_m31.sh <runname> [steps] [patience] [layer]
RUN=${1:-runA}; STEPS=${2:-4608}; PAT=${3:-2}; LAYER=${4:--1}
O=$HOME/coding-rna/outputs/M32-C1-PLANTCAD2; mkdir -p $O/$RUN
echo "{\"event\":\"queued\",\"time\":\"$(date -Is)\"}" >> $O/$RUN/queue.jsonl
gone=0
while [ $gone -lt 3 ]; do        # M31 processes absent in 3 consecutive checks, 5 min apart (bridges the short gaps between M31 stage steps)
  if pgrep -f 'M31-C1-NT[V]2' > /dev/null; then gone=0; else gone=$((gone+1)); fi
  [ $gone -lt 3 ] && sleep 300
done
echo "{\"event\":\"m31_gone_start_m32\",\"time\":\"$(date -Is)\"}" >> $O/$RUN/queue.jsonl
cd $HOME/coding-rna-reorg/experiments/M32-C1-PLANTCAD2/scripts
exec bash c1_stage.sh $RUN $STEPS $PAT $LAYER
