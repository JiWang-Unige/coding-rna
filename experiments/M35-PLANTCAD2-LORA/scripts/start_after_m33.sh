#!/bin/bash
# Queue driver (detached on ont): GPU rule for M35 = start only when ~/coding-rna/outputs/M33-NTV2-FOLLOWUP/M33_DONE exists AND no process whose command line contains M33-NTV2-FOLLOWUP
# is left, in 2 consecutive checks 5 min apart.  Then run the GPU smoke (smoke.sh) and STOP (the full run is launched after the smoke has been reviewed).
O=$HOME/coding-rna/outputs/M35-PLANTCAD2-LORA; mkdir -p $O
echo "{\"event\":\"queued\",\"time\":\"$(date -Is)\"}" >> $O/queue.jsonl
ok=0
while [ $ok -lt 2 ]; do
  if [ -f $HOME/coding-rna/outputs/M33-NTV2-FOLLOWUP/M33_DONE ] && ! pgrep -f 'M33-NTV2-FOLLOW[U]P' > /dev/null; then ok=$((ok+1)); else ok=0; fi
  [ $ok -lt 2 ] && sleep 300
done
echo "{\"event\":\"m33_gate_open_start_smoke\",\"time\":\"$(date -Is)\"}" >> $O/queue.jsonl
cd $HOME/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA/scripts
bash smoke.sh >> $O/queue_smoke.out 2>&1
echo "{\"event\":\"smoke_finished\",\"time\":\"$(date -Is)\"}" >> $O/queue.jsonl
