#!/bin/bash
# PANEL-P6 GPU queue: SMOKE -> ANNEVO -> Tiberius -> Helixer -> Carbon-A -> Vipsania (each tool over the core-10 genomes, smallest first).
# SMOKE: every tool of P6_SMOKE_TOOLS end to end on a 2-Mb Arabidopsis TRAINING slice (not scored); any failure stops the queue before the production runs
# ($P6/SMOKE_OK when passed).  Every GPU step waits at the gate of p6_lib.sh ($P6/GO, no TE GPU work, no other coding-rna GPU job, enough free memory and disk).
# Single instance; resumable (finished genomes are skipped).  P6_STEPS / P6_SMOKE_TOOLS select a subset (e.g. win3090: the three native tools only, p6_win3090.sh).
# Status: tail -n 5 ~/coding-rna/outputs/PANEL-P6/events.jsonl      Stop between steps: touch ~/coding-rna/outputs/PANEL-P6/STOP
# Launch: cd ~/coding-rna-reorg/experiments/PANEL-P6/scripts && (setsid nohup bash p6_queue.sh > ~/coding-rna/logs/p6-queue.log 2>&1 < /dev/null &)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HERE/p6_lib.sh
exec 9>$P6/.queue.lock
flock -n 9 || { echo "p6_queue already running"; exit 3; }
STEPS=${P6_STEPS:-SMOKE ANNEVO Tiberius Helixer CarbonA Vipsania}
SMOKE_TOOLS=${P6_SMOKE_TOOLS:-ANNEVO Tiberius Helixer CarbonA Vipsania}
ev event=queue_start pid=$$ host=$(hostname) "steps=$STEPS" "smoke_tools=$SMOKE_TOOLS"
for s in $STEPS; do
  case $s in
    SMOKE)
      if [ ! -e $P6/SMOKE_OK ]; then
        ok=1
        for t in $SMOKE_TOOLS; do
          case $t in
            ANNEVO|Tiberius|Helixer) P6_SMOKE=1 bash $HERE/p6_native.sh $t;;
            CarbonA) P6_SMOKE=1 bash $HERE/p6_carbona.sh;;
            Vipsania) P6_SMOKE=1 bash $HERE/p6_vipsania.sh;;
          esac || { ok=0; ev event=smoke_failed tool=$t; break; }
        done
        [ $ok = 1 ] && { echo "$(date -Is) $SMOKE_TOOLS" > $P6/SMOKE_OK; ev event=smoke_ok "tools=$SMOKE_TOOLS"; }
        [ $ok = 1 ]
      fi;;
    ANNEVO|Tiberius|Helixer) bash $HERE/p6_native.sh $s;;
    CarbonA) bash $HERE/p6_carbona.sh;;
    Vipsania) bash $HERE/p6_vipsania.sh;;
  esac
  rc=$?
  [ $rc -eq 0 ] || { ev event=queue_aborted step=$s rc=$rc; exit $rc; }
done
date -Is > $P6/P6_GPU_DONE_$(hostname)
ev event=queue_done
exit 0
