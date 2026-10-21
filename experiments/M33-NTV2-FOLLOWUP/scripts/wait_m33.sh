#!/bin/bash
# Event-driven waiter for the M33 stages (Mac, background).  Returns on: gpu_inference.done | inference+postq both gone | new error lines | status change in
# interruptions | postq 'fail' | heartbeat.  usage: wait_m33.sh [heartbeat_s=3000] [interval_s=1200] [stop_when: gpu|postq|any]
HB=${1:-3000}; IV=${2:-1200}; MODE=${3:-gpu}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
P="bash ~/coding-rna-reorg/experiments/M33-NTV2-FOLLOWUP/scripts/probe_m33.sh"
t0=$(date +%s); e0=""
while true; do
  sleep $IV
  st=$(Q "$P"); [ -z "$st" ] && continue
  l1=$(echo "$st" | head -1); ni=${l1%% *}; rest=${l1#* }; np=${rest%% *}; rest=${rest#* }; g=${rest%% *}; ne=${rest#* }
  [ -z "$e0" ] && e0=$ne
  el=$(( $(date +%s) - t0 )); reason=""
  [ "$ne" != "$e0" ] && reason="new_error_lines"
  echo "$st" | grep -q '"fail"' && reason="postq_fail"
  [ "$MODE" = "gpu" ] && [ "$g" = "1" ] && reason="gpu_inference_done"
  [ "$MODE" = "postq" ] && [ "$np" = "0" ] && reason="postq_exited"
  [ "$ni" = "0" ] && [ "$g" != "1" ] && reason="inference_exited_without_done"
  [ -z "$reason" ] && [ $el -ge $HB ] && reason="heartbeat_${el}s"
  [ -n "$reason" ] && { echo "EVENT: $reason"; echo "$st"; exit 0; }
done
