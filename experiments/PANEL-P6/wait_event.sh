#!/bin/bash
# Event-driven waiter (runs on the Mac, background).  Returns when the P6 queue logs an event other than gate_wait, the queue process is gone, a given ont file appears,
# or the heartbeat expires (default 50 min; checks every 10 min).   usage: wait_event.sh [WATCH_FILE_ON_ONT] [heartbeat_s=3000] [interval_s=600]
W=${1:-/nonexistent}; HB=${2:-3000}; IV=${3:-600}; t0=$(date +%s)
Q() { ssh baobab "ssh lab \"$1\" < /dev/null" 2>/dev/null; }
W3() { ssh -o ConnectTimeout=30 win3090-wsl "$1" < /dev/null 2>/dev/null; }      # win3090 queue (native tools)
n0=$(Q "grep -vc gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl")
w0=$(W3 "grep -vc gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl")
while true; do
  sleep $IV
  out=$(Q "grep -vc gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl; pgrep -f 'bash p6_queue.sh' | wc -l; ls $W 2>/dev/null | wc -l")
  n=$(echo "$out" | sed -n 1p); alive=$(echo "$out" | sed -n 2p); seen=$(echo "$out" | sed -n 3p)
  el=$(( $(date +%s) - t0 )); reason=""
  w=$(W3 "grep -vc gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl")
  [ -n "$w" ] && [ -n "$w0" ] && [ "$w" != "$w0" ] && reason="win3090_event"
  [ -n "$n" ] && [ "$n" != "$n0" ] && reason="new_event"
  [ "$alive" = "0" ] && reason="queue_exited"
  [ "$seen" != "0" ] && [ -n "$seen" ] && reason="watch_file:$W"
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason"
    Q "grep -v gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl | tail -n 4 | cut -c1-400; tail -n 1 ~/coding-rna/outputs/PANEL-P6/events.jsonl | cut -c1-600; nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader; df -h ~ | tail -1"
    W3 "grep -v gate_wait ~/coding-rna/outputs/PANEL-P6/events.jsonl | tail -n 3 | cut -c1-300; pgrep -fc 'bash /home/ji/coding-rna-reorg/experiments/PANEL-P6/scripts/p6_[q]ueue.sh'" | sed 's/^/win3090: /'
    exit 0
  fi
done
