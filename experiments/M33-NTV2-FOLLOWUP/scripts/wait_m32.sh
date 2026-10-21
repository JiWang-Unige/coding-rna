#!/bin/bash
# Event-driven waiter (runs on the Mac, background): returns when ANY of:
#   M32_DONE exists AND no M32 process | M32 processes vanished without M32_DONE (after having been seen) | my own M33 job (if any) exits / errors | 50-min heartbeat.
# Checks over ssh every IV seconds (default 20 min).  usage: wait_m32.sh [heartbeat_s=3000] [interval_s=1200]
HB=${1:-3000}; IV=${2:-1200}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
t0=$(date +%s); seen=0
while true; do
  st=$(Q "bash ~/coding-rna-reorg/experiments/M33-NTV2-FOLLOWUP/scripts/probe_m32.sh")
  if [ -z "$st" ]; then sleep 300; el=$(( $(date +%s) - t0 )); [ $el -ge $HB ] && { echo "EVENT: heartbeat_probe_unreachable_${el}s"; exit 0; }; continue; fi   # ssh/network failure: not evidence about M32
  n=$(echo "$st" | head -1 | cut -d' ' -f1); d=$(echo "$st" | head -1 | cut -d' ' -f2)
  [ "${n:-0}" -gt 0 ] && seen=1
  el=$(( $(date +%s) - t0 )); reason=""
  [ "$d" = "1" ] && [ "${n:-0}" = "0" ] && reason="M32_DONE_and_no_M32_process"
  [ "$d" = "1" ] && [ "${n:-0}" != "0" ] && reason="M32_DONE_but_M32_processes_still_listed"
  [ -z "$reason" ] && [ $el -ge $HB ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then echo "EVENT: $reason"; echo "$st"; exit 0; fi
  sleep $IV
done
