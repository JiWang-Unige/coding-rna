#!/bin/bash
# Event-driven waiter (runs on the Mac): returns when the log shows ALL_DONE / a MISMATCH / curl error, the process is gone, or a 50-min heartbeat.
# usage: wait_event.sh <log_prefix_on_ont> [heartbeat_seconds] [check_interval_seconds]
LOG=$1; HB=${2:-3000}; IV=${3:-60}; t0=$(date +%s)
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
while true; do
  sleep $IV
  out=$(Q "tail -n 3 $LOG.out; echo ERR; tail -n 3 $LOG.err; echo ALIVE; pgrep -f m34_download | wc -l")
  el=$(( $(date +%s) - t0 )); reason=""
  echo "$out" | grep -q ALL_DONE && reason="all_done"
  echo "$out" | grep -q MISMATCH && reason="md5_mismatch"
  [ -n "$(echo "$out" | sed -n '/^ERR/,/^ALIVE/p' | grep -v -e '^ERR' -e '^ALIVE')" ] && reason="stderr_output"
  [ "$(echo "$out" | tail -1)" = "0" ] && [ -z "$reason" ] && reason="process_exited"
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then echo "EVENT: $reason"; echo "$out"; Q "du -sh ~/coding-rna/data/m34_train/raw"; exit 0; fi
done
