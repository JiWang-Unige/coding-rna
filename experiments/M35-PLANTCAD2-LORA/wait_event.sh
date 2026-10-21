#!/bin/bash
# Event-driven waiter (runs on the Mac, start with run_in_background): returns when ANY of: sentinel (SMOKE_DONE / M35_DONE / STOP_REQUEST) appears | stage/metrics.jsonl gains a line | >=3 new
# interruptions | Traceback/Killed in an M35 err log | all M35 processes gone after having started | heartbeat.  Checks every IV seconds (default 20 min).
# usage: wait_event.sh <run (e.g. runB)> [heartbeat_seconds=3000] [check_interval_seconds=1200]
RUN=${1:-runB}; HB=${2:-3000}; IV=${3:-1200}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
O='/home/ji/coding-rna/outputs/M35-PLANTCAD2-LORA'; R=$O/$RUN
n0=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i0=$(Q "cat $O/smoke/interruptions.jsonl $R/interruptions.jsonl 2>/dev/null | grep -c interruption")
d0=$(Q "ls $O/smoke/SMOKE_DONE $O/M35_DONE 2>/dev/null | wc -l")
TB() { Q "cat ~/coding-rna/logs/m35-*.err ~/coding-rna/logs/c1-infer-*.err 2>/dev/null | grep -c -E 'Traceback|Killed'"; }
tb0=$(TB)
t0=$(date +%s); seen=0
while true; do
  sleep $IV
  alive=$(Q "pgrep -f 'M35-PLANTCAD2-LOR[A]' | wc -l")
  [ "${alive:-0}" != "0" ] && seen=1
  n=$(Q "wc -l < $R/metrics.jsonl 2>/dev/null || echo 0"); i=$(Q "cat $O/smoke/interruptions.jsonl $R/interruptions.jsonl 2>/dev/null | grep -c interruption")
  d=$(Q "ls $O/smoke/SMOKE_DONE $O/M35_DONE 2>/dev/null | wc -l")
  tb=$(TB)
  el=$(( $(date +%s) - t0 )); reason=""
  [ "${d:-0}" != "${d0:-0}" ] && reason="sentinel_appeared"
  [ "${n:-0}" != "${n0:-0}" ] && reason="new_metrics_line"
  [ $(( ${i:-0} - ${i0:-0} )) -ge 3 ] && reason="3+_new_interruptions"
  [ "${tb:-0}" -gt "${tb0:-0}" ] && reason="new_traceback_or_killed_in_err_logs"
  [ "$seen" = "1" ] && [ "${alive:-0}" = "0" ] && reason="m35_processes_exited"
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason"
    Q "tail -n 3 $O/queue.jsonl 2>/dev/null; tail -n 4 $O/smoke/smoke.jsonl 2>/dev/null | cut -c1-200; echo ---; tail -n 2 $R/metrics.jsonl 2>/dev/null | cut -c1-400; echo ---; tail -n 3 $R/interruptions.jsonl 2>/dev/null | cut -c1-300; tail -n 2 $R/stage_log.jsonl 2>/dev/null; echo ---; tail -n 2 $R/training.jsonl 2>/dev/null | cut -c1-460; tail -n 2 $O/smoke/training.jsonl 2>/dev/null | cut -c1-460; ls $O/M33*/ 2>/dev/null | head -0; ls /home/ji/coding-rna/outputs/M33-NTV2-FOLLOWUP/M33_DONE 2>&1 | tail -1; nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader"
    exit 0
  fi
done
