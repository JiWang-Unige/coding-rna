#!/bin/bash
# Keeps the k50F -> k25U chain (m34_after_m36.sh) alive across GPU-memory waits.  m34_supervise.sh gives up after 80 waits of 120 s (< 7.5 GB free) and m34_stage.sh then
# exits 9; every stage is resumable, so when the chain process is gone without all_done, this waits until the GPU has >= 7.5 GB free for 3 consecutive checks
# (120 s apart) and starts m34_after_m36.sh again.  At most 6 relaunches; events go to after_m36.jsonl.
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
LOG=$O/after_m36.jsonl
exec 8>$O/.relaunch_after_m36.lock
flock -n 8 || { echo "relauncher already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ev relauncher_started
n=0
while [ $n -lt 6 ]; do
  while pgrep -u $(id -un) -f 'm34_after_m3[6]' > /dev/null; do sleep 120; done
  grep -q '"event":"all_done"' $LOG && { ev relauncher_exit ",\"reason\":\"all_done\""; exit 0; }
  ok=0
  while [ $ok -lt 3 ]; do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "$free" -ge 7500 ]; then ok=$((ok+1)); else ok=0; fi
    sleep 120
  done
  n=$((n+1))
  ev relaunch ",\"n\":$n,\"free_mib\":$free"
  setsid nohup bash m34_after_m36.sh >> $HOME/coding-rna/logs/m34-after-m36.out 2>&1 < /dev/null &
  sleep 60
done
ev relauncher_exit ",\"reason\":\"relaunch budget exhausted\""
