#!/bin/bash
# Event-driven waiter for the M34 queue (runs on the Mac, background).  Returns when ANY of: a new queue.jsonl / stage_log.jsonl / metrics.jsonl line,
# M34_STAGE1_DONE, a FATAL / failed event, 3+ new interruptions, a Traceback in the newest train log, the queue process gone (without STAGE1_DONE), or the heartbeat
# (default 50 min).  Checks every IV seconds (default 20 min).  usage: wait_stage.sh [heartbeat_s=3000] [interval_s=1200]
HB=${1:-3000}; IV=${2:-1200}
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
O='~/coding-rna/outputs/M34-NTV2-SCALEUP'
snap() { Q "echo \$(cat $O/queue.jsonl 2>/dev/null | wc -l) \$(cat $O/*/stage_log.jsonl 2>/dev/null | wc -l) \$(cat $O/*/metrics.jsonl 2>/dev/null | wc -l) \$(cat $O/*/interruptions.jsonl 2>/dev/null | grep -c '\"interruption\"') \$(ls $O/M34_STAGE1_DONE 2>/dev/null | wc -l) \$(pgrep -f 'm34_queue.sh' | wc -l) \$(cat $O/*/interruptions.jsonl 2>/dev/null | grep -c FATAL)"; }
read -r q0 s0 m0 i0 d0 p0 f0 <<< "$(snap)"
t0=$(date +%s)
while true; do
  sleep $IV
  read -r q s m i d p f <<< "$(snap)"
  el=$(( $(date +%s) - t0 )); reason=""
  [ -z "$q" ] && continue
  [ "$d" != "0" ] && reason="STAGE1_DONE"
  [ "${f:-0}" -gt "${f0:-0}" ] && reason="FATAL_event"
  [ "$q" != "$q0" ] && reason="new_queue_event"
  [ "$s" != "$s0" ] && reason="new_stage_event"
  [ "$m" != "$m0" ] && reason="new_metrics_line"
  [ $(( ${i:-0} - ${i0:-0} )) -ge 3 ] && reason="3+_new_interruptions"
  [ "$p" = "0" ] && [ "$d" = "0" ] && reason="queue_process_gone"
  [ -z "$reason" ] && [ $el -ge $HB ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason"
    Q "tail -n 3 $O/queue.jsonl 2>/dev/null; echo ---; for a in smoke T2 T10 pc2_smoke pc2_T2 pc2_T10; do [ -d $O/\$a ] && { echo \$a; tail -n 2 $O/\$a/stage_log.jsonl 2>/dev/null | cut -c1-300; tail -n 1 $O/\$a/metrics.jsonl 2>/dev/null | cut -c1-260; tail -n 2 $O/\$a/interruptions.jsonl 2>/dev/null | cut -c1-300; tail -n 1 $O/\$a/training.jsonl 2>/dev/null | cut -c1-260; }; done; nvidia-smi --query-gpu=memory.used --format=csv,noheader; bash ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts/m34_gate_status.sh 2>/dev/null | tail -2"
    exit 0
  fi
done
