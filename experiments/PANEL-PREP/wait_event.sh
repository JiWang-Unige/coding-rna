#!/bin/bash
# Event-driven waiter (runs on the Mac, background): returns when the log shows the DONE marker, a MISMATCH/Traceback/Error, the process pattern is gone,
# or the 50-min heartbeat.  usage: wait_event.sh <log_on_ont> <pgrep_pattern> [done_marker=ALL_DONE] [heartbeat_s=3000] [interval_s=120]
LOG=$1; PAT=$2; DONE=${3:-ALL_DONE}; HB=${4:-3000}; IV=${5:-120}; t0=$(date +%s)
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
while true; do
  sleep $IV
  out=$(Q "tail -n 4 $LOG; echo ALIVE; pgrep -f '$PAT' | wc -l")
  el=$(( $(date +%s) - t0 )); reason=""
  echo "$out" | grep -q "$DONE" && reason="done"
  echo "$out" | grep -q -e MISMATCH -e Traceback -e "Error" -e Killed -e FAILED && reason="error_marker"
  [ "$(echo "$out" | tail -1)" = "0" ] && [ -z "$reason" ] && reason="process_exited"
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then echo "EVENT: $reason"; echo "$out"; Q "du -sh ~/coding-rna/data/panel 2>/dev/null; df -h ~ | tail -1"; exit 0; fi
done
