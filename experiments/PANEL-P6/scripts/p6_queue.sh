#!/bin/bash
# PANEL-P6 GPU queue on ont: ANNEVO -> Tiberius -> Helixer -> Carbon-A -> Vipsania (each tool over the core-10 genomes, smallest first).  Every GPU step waits at the gate
# of p6_lib.sh ($P6/GO, no TE GPU process, no other coding-rna GPU job, enough free memory and disk).  Single instance; resumable (finished genomes are skipped).
# Status: tail -n 5 ~/coding-rna/outputs/PANEL-P6/events.jsonl      Stop between steps: touch ~/coding-rna/outputs/PANEL-P6/STOP
# Launch: cd ~/coding-rna-reorg/experiments/PANEL-P6/scripts && (setsid nohup bash p6_queue.sh > ~/coding-rna/logs/p6-queue.log 2>&1 < /dev/null &)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HERE/p6_lib.sh
exec 9>$P6/.queue.lock
flock -n 9 || { echo "p6_queue already running"; exit 3; }
ev event=queue_start pid=$$ "steps=${P6_STEPS:-ANNEVO Tiberius Helixer CarbonA Vipsania}"
for s in ${P6_STEPS:-ANNEVO Tiberius Helixer CarbonA Vipsania}; do
  case $s in
    ANNEVO|Tiberius|Helixer) bash $HERE/p6_native.sh $s;;
    CarbonA) bash $HERE/p6_carbona.sh;;
    Vipsania) bash $HERE/p6_vipsania.sh;;
  esac
  rc=$?
  [ $rc -eq 0 ] || { ev event=queue_aborted step=$s rc=$rc; exit $rc; }
done
date -Is > $P6/P6_GPU_DONE
ev event=queue_done
exit 0
