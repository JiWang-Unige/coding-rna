#!/bin/bash
# Event-driven waiter (runs on the Mac): returns when ANY of: stage process gone | new metrics.jsonl line | new interruption event | Traceback in err logs | 50-min heartbeat.
# usage: wait_event.sh <run> [heartbeat_seconds] [check_interval_seconds]
RUN=${1:-runA}; HB=${2:-3000}; IV=${3:-1200}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
R='~/coding-rna/outputs/M31-C1-NTV2/'$RUN
n0=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i0=$(Q "grep -c interruption $R/interruptions.jsonl 2>/dev/null" | head -1); i0=${i0:-0}
t0=$(date +%s)
while true; do
  sleep $IV
  alive=$(Q "pgrep -f c1_stage | wc -l")
  n=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i=$(Q "grep -c interruption $R/interruptions.jsonl 2>/dev/null" | head -1); i=${i:-0}
  tb=$(Q "grep -l -E 'Traceback|Killed' ~/coding-rna/logs/m31-train-$RUN-*.err 2>/dev/null | head -1")
  el=$(( $(date +%s) - t0 ))
  reason=""
  [ "${alive:-0}" = "0" ] && reason="stage_process_exited"
  [ "${n:-0}" != "${n0:-0}" ] && reason="new_metrics_line"
  [ $(( ${i:-0} - ${i0:-0} )) -ge 3 ] && reason="3+_new_interruptions"
    [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason"
    Q "tail -n 2 $R/metrics.jsonl 2>/dev/null; echo ---; tail -n 3 $R/interruptions.jsonl; echo ---; tail -n 2 $R/stage_log.jsonl 2>/dev/null; echo ---; tail -n 2 $R/training.jsonl 2>/dev/null | cut -c1-420; nvidia-smi --query-gpu=memory.used --format=csv,noheader"
    exit 0
  fi
done
