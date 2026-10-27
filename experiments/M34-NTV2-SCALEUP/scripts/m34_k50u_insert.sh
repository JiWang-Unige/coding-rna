#!/bin/bash
# Insert M34 arm T10_k50U (pc2) into the GPU gap while the M36 queue waits for its CPU seeds (user decision gpu-gap, 2026-10-27).
# Start condition (checked every 120 s): M36 arm S2_eschscholzia_californica has finished (DONE_512.json), the M36 S1 seeds of the next species are NOT ready
# (READY_dryas_octopetala missing), and no coding-rna GPU process for 3 consecutive checks.  If the seeds become ready first, nothing is inserted (k50U keeps its
# place after M36) and the script exits.  While k50U runs, the M36 queue process is frozen with SIGSTOP (its CPU seed queue is a separate process and keeps going) and
# resumed with SIGCONT when k50U ends or this script exits for any reason -- so the two coding-rna GPU jobs never overlap (k50U's validation decode is CPU-only, a GPU
# gate alone would let M36 start in between).
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
M36=$HOME/coding-rna/outputs/M36-BOOTSTRAP-ADAPT
LOG=$O/k50u_insert.jsonl
exec 7>$O/.k50u_insert.lock
flock -n 7 || { echo "k50u insert already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
source $HOME/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts/m36_gate.sh      # m36_gpu_blockers (coding-rna GPU pids of this user)
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $LOG; }
ev queued
ok=0; n=0
while :; do
  if [ -f $M36/seeds/READY_dryas_octopetala ]; then ev skipped_seeds_ready; exit 0; fi
  if [ -f $M36/train/S2_eschscholzia_californica/DONE_512.json ] && [ -z "$(m36_gpu_blockers)" ]; then ok=$((ok+1)); else ok=0; fi
  [ $ok -ge 3 ] && break
  [ $((n % 30)) -eq 0 ] && ev waiting ",\"ok\":$ok"
  n=$((n+1)); sleep 120
done
M36PID=$(pgrep -u $(id -un) -f 'm36_queue.s[h]' | head -1)
[ -n "$M36PID" ] || { ev no_m36_queue_process; }
resume() { [ -n "$M36PID" ] && kill -CONT $M36PID 2>/dev/null && ev m36_resumed ",\"pid\":$M36PID"; }
trap resume EXIT
if [ -n "$M36PID" ]; then kill -STOP $M36PID && ev m36_paused ",\"pid\":$M36PID"; fi
if [ -f $M36/seeds/READY_dryas_octopetala ]; then ev skipped_seeds_ready_late; exit 0; fi
export M34_BACKEND=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}') M34_LORA=on
source ./m34_env.sh
ev k50u_start ",\"backend\":\"$M34_BACKEND\",\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
bash m34_stage.sh T10_k50U; rc=$?
ev k50u_end ",\"rc\":$rc"
exit $rc
