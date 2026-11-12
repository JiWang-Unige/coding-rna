#!/bin/bash
# One-off GPU smoke of the PlantCAD2 (pc2) backend on the M39 insect data (decision pc2-insect-transfer, user 2026-11-11: smoke first, ~1 GPU h, in a gap
# that TE leaves on the ont GPU).  = `M39_BACKEND=pc2 m39_smoke.sh` (I10_T 10 + resume 10 steps, k50U 3, k50F 3, 300-kb slice inference + scoring)
# -> ~/coding-rna/outputs/M39-INSECT-C2/pc2_smoke/{summary.json, SMOKE_OK | SMOKE_FAIL}.  Touches no NT-v2 output (smoke/, I10_*/).
#   start condition: no TE_GPU_CLAIM, no TE process on the GPU (te_on_gpu() of m39_after_k25u.sh), and the GPU has NO compute process at all and >= 7.5 GB free,
#                    for 6 consecutive checks 5 min apart (30-min TE quiet period, as the P6 gate since 11-11).  It never starts while the P6 queue or the M34 queue has something on the GPU.
#   yield to TE: if TE_GPU_CLAIM appears or a TE process shows up on the GPU while the smoke runs, the smoke's own process group is killed (event
#                yield_to_te) and the gate waits again (the rerun starts the smoke from scratch).  Nothing else is ever killed.
#   gives up after M39_GATE_HOURS (default 72) h of waiting (event gave_up).
# Events: ~/coding-rna/outputs/M39-INSECT-C2/pc2_smoke/gate.jsonl ; nvidia-smi before the start: pc2_smoke/nvidia_smi_start_<try>.txt
# Launch (ont): cd ~/coding-rna-reorg/experiments/M39-INSECT-C2/scripts && (setsid nohup bash m39_pc2_smoke_gate.sh > ~/coding-rna/logs/m39-pc2-smoke-gate.out 2>&1 < /dev/null &)
# Status:       tail -n 3 ~/coding-rna/outputs/M39-INSECT-C2/pc2_smoke/gate.jsonl; cat ~/coding-rna/outputs/M39-INSECT-C2/pc2_smoke/summary.json
cd $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts
export M39_BACKEND=${M39_BACKEND:-pc2}          # also used for the one-hot floor (M39_BACKEND=onehot -> oh_smoke/; protocol 32 arm L)
source ./m39_env.sh
O=$OUTROOT/${ARMPFX}smoke
mkdir -p $O
exec 7>$O/.gate.lock
flock -n 7 || { echo "pc2 smoke gate already running"; exit 3; }
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $O/gate.jsonl; }
te_on_gpu() {   # identical to m39_after_k25u.sh / m34_after_p6.sh
  local p
  [ -e $HOME/ab-initio-TE/outputs/TE_GPU_CLAIM ] && return 0
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -qE 'ab-initio-TE|te_benchmark|TE-GLM|/work/te' && return 0
    readlink /proc/$p/cwd 2>/dev/null | grep -qE 'ab-initio-TE|TE-GLM|/work/te' && return 0
    ps -o comm= -p $(ps -o ppid= -p $p 2>/dev/null) 2>/dev/null | grep -qiE 'apptainer|singularity' && return 0
  done
  return 1
}
gpu_idle() {    # no compute process at all and >= 7.5 GB free.  M39_GATE_SHARE=1: other (non-TE) processes allowed, only >= 7.5 GB free
  local free napps         # (share mode: maintenance session 11-11 OK'd sharing the GPU with P6 / k50U_s2; TE still blocks via te_on_gpu)
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  napps=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -c .)
  { [ "${M39_GATE_SHARE:-0}" = 1 ] || [ "${napps:-1}" -eq 0 ]; } && [ "${free:-0}" -ge 7500 ]
}
DEADLINE=$(( $(date +%s) + ${M39_GATE_HOURS:-72} * 3600 ))
ev gate_start ",\"backend\":\"$M39_BACKEND\",\"share\":${M39_GATE_SHARE:-0}"
for try in 1 2 3 4 5 6; do
  [ -f $O/SMOKE_OK ] && { ev smoke_ok_exists; exit 0; }
  ok=0; n=0
  while [ $ok -lt 6 ]; do
    if ! te_on_gpu && gpu_idle; then ok=$((ok+1)); else ok=0; fi
    [ $((n % 12)) -eq 0 ] && [ $ok -eq 0 ] && ev waiting ",\"te\":$(te_on_gpu && echo 1 || echo 0),\"gpu_apps\":$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -c .)"
    n=$((n+1))
    [ $(date +%s) -gt $DEADLINE ] && { ev gave_up; exit 1; }
    [ $ok -lt 6 ] && sleep 300
  done
  nvidia-smi > $O/nvidia_smi_start_$try.txt 2>&1
  ev smoke_start ",\"try\":$try,\"gpu\":\"$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader | tr -d ' ')\""
  setsid bash m39_smoke.sh >> $HOME/coding-rna/logs/m39-pc2-smoke.out 2>&1 < /dev/null &
  SPID=$!
  yielded=0
  while kill -0 $SPID 2>/dev/null; do
    if te_on_gpu; then
      ev yield_to_te ",\"try\":$try"
      kill -TERM -- -$SPID 2>/dev/null; sleep 20; kill -KILL -- -$SPID 2>/dev/null
      yielded=1; break
    fi
    sleep 30
  done
  wait $SPID 2>/dev/null; rc=$?
  nvidia-smi > $O/nvidia_smi_end_$try.txt 2>&1
  if [ $yielded -eq 1 ]; then rm -rf $O/t $O/k50U $O/k50F $O/infer_slice $O/metrics.jsonl; continue; fi    # clean rerun (the smoke checks steps 1..20 exactly once)
  if [ -f $O/SMOKE_OK ]; then ev smoke_ok ",\"try\":$try"; exit 0; fi
  ev smoke_failed ",\"try\":$try,\"rc\":$rc"; exit 9          # a real failure is not retried: inspect pc2_smoke/smoke_log.jsonl and the logs
done
ev gave_up ",\"reason\":\"6 tries (TE reclaimed the GPU each time)\""
exit 1
