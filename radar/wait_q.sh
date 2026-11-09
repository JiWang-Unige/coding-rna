#!/bin/bash
# Event-driven waiter for the ont GPU queue (run on the Mac in background by the dashboard-maintainer session).  Returns when ANY of:
# new queue / panel / M38 event (counted by ~/coding-rna/tools/q_status.sh) | a sentinel appears | a stage metrics/stage_log line is added |
# Traceback/Killed/FATAL in the watched err logs or a failed posterior task | queue process gone after having been seen |
# free disk on / below DISK GB (arg 3, 0 = off) | heartbeat (default 50 min).  Checks every IV s.
# usage: bash radar/wait_q.sh [heartbeat_s=3000] [check_s=900] [disk_gb=0]
HB=${1:-3000}; IV=${2:-900}; DISK=${3:-0}
CM="-o ControlMaster=auto -o ControlPath=$HOME/.ssh/cm-waitq-%r@%h-%p -o ControlPersist=20m"   # reuse one baobab connection (11-09: the login node reset frequent new connections)
Q() { ssh $CM -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ServerAliveCountMax=4 baobab "ssh -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ServerAliveCountMax=4 lab \"$1\" < /dev/null" 2>/dev/null; }
S() { Q "bash ~/coding-rna/tools/q_status.sh counts" | tr '\n' ' '; }
for i in 1 2 3; do read q0 s0 m0 e0 a0 <<<"$(S)"; [ -n "$a0" ] && break; sleep 60; done; t0=$(date +%s); seen=0
while true; do
  sleep $IV
  read q s m e a <<<"$(S)"
  el=$(( $(date +%s) - t0 )); reason=""
  if [ -z "$a" ]; then
    miss=$(( ${miss:-0} + 1 ))
    if [ $miss -ge 3 ]; then echo "EVENT: ssh_unreachable_3x  ($(date '+%F %T'))"; exit 0; fi
    if [ $el -ge $HB ]; then echo "EVENT: heartbeat_${el}s_but_ssh_failed  ($(date '+%F %T'))"; exit 0; fi
    continue
  fi
  miss=0
  [ "${a:-0}" != "0" ] && seen=1
  [ "${q:-0}" != "${q0:-0}" ] && reason="new_queue_or_panel_event"
  [ "${s:-0}" != "${s0:-0}" ] && reason="sentinel"
  [ "${m:-0}" != "${m0:-0}" ] && reason="new_stage_metrics_or_log_line"
  [ "${e:-0}" -gt "${e0:-0}" ] && reason="new_error"
  [ "$seen" = "1" ] && [ "${a:-0}" = "0" ] && reason="queue_process_exited"
  if [ "$DISK" -gt 0 ] && [ -z "$reason" ]; then
    fr=$(Q "df -BG --output=avail / | tail -1" | tr -dc 0-9)
    [ -n "$fr" ] && [ "$fr" -lt "$DISK" ] && reason="disk_free_${fr}G_below_${DISK}G"
  fi
  [ $el -ge $HB ] && [ -z "$reason" ] && reason="heartbeat_${el}s"
  if [ -n "$reason" ]; then
    echo "EVENT: $reason  ($(date '+%F %T'))"
    Q "bash ~/coding-rna/tools/q_status.sh details"
    exit 0
  fi
done
