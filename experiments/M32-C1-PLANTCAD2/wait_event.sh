#!/bin/bash
# Event-driven waiter (runs on the Mac): returns when ANY of: queue/stage process gone (after it started) | new metrics.jsonl line | 3+ new interruptions | Traceback in err logs | heartbeat.
# usage: wait_event.sh <run> [heartbeat_seconds=3000] [check_interval_seconds=1200]
RUN=${1:-runA}; HB=${2:-3000}; IV=${3:-1200}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
R='~/coding-rna/outputs/M32-C1-PLANTCAD2/'$RUN
n0=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i0=$(Q "grep -c interruption $R/interruptions.jsonl 2>/dev/null || echo 0")
t0=$(date +%s)
while true; do
  sleep $IV
  alive=$(Q "pgrep -f 'start_after_m3[1]|M32-C1-PLANTCAD[2]' | wc -l")
  n=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i=$(Q "grep -c interruption $R/interruptions.jsonl 2>/dev/null || echo 0")
  tb=$(Q "grep -l -E 'Traceback|Killed' ~/coding-rna/logs/m32-train-$RUN-*.err 2>/dev/null | head -1")
  el=$(( $(date +%s) - t0 ))
  reason=""
  [ "${alive:-0}" = "0" ] && reason="queue_and_stage_processes_exited"
  [ "${n:-0}" != "${n0:-0}" ] && reason="new_metrics_line"
  [ $(( ${i:-0} - ${i0:-0} )) -ge 3 ] && reason="3+_new_interruptions"
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason"
    Q "tail -n 2 $R/metrics.jsonl 2>/dev/null | cut -c1-400; echo ---; tail -n 3 $R/queue.jsonl 2>/dev/null; tail -n 3 $R/interruptions.jsonl 2>/dev/null | cut -c1-300; echo ---; tail -n 2 $R/stage_log.jsonl 2>/dev/null; echo ---; tail -n 2 $R/training.jsonl 2>/dev/null | cut -c1-420; nvidia-smi --query-gpu=memory.used --format=csv,noheader"
    exit 0
  fi
done
