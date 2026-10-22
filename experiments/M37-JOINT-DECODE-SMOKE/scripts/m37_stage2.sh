#!/bin/bash
# Stage 2 (after stage 1): finish B inference if the 20-min watchdog cut it (resume is per strand AND per window batch via the encode cache
# of m37_infer.py; each run again <= 20 min GPU; up to 12 attempts, 2026-10-21 23:17 protocol §8),
# then joint decoding (lambda 0 / 0.25 / 0.5), unit tests (incl. lambda=0 identity on the real chromosome), analysis.  Marker: $O/stage2.done
O=$HOME/coding-rna/outputs/M37-JOINT-DECODE-SMOKE; L=$HOME/coding-rna/logs
S=$HOME/coding-rna-reorg/experiments/M37-JOINT-DECODE-SMOKE; PY=$HOME/envs/generanno/bin/python
cd $S/scripts
while [ ! -f $O/stage1.done ]; do sleep 60; done
for attempt in $(seq 3 14); do
  ls $O/inferB/DONE_*.json >/dev/null 2>&1 && break
  pgrep -f m37_infer.py >/dev/null && { sleep 60; continue; }
  echo "{\"event\":\"resume_attempt\",\"attempt\":$attempt,\"time\":\"$(date -Is)\"}" >> $O/inferB/gpu_log.jsonl
  bash m37_gpu_infer.sh; rc=$?
  [ $rc -eq 2 ] && sleep 300          # GPU memory below 8 GB: wait and retry
done
ls $O/inferB/DONE_*.json >/dev/null 2>&1 || { echo "inferB_not_done $(date -Is)" >> $O/stage2.log; exit 1; }
echo "inferB done $(date -Is)" >> $O/stage2.log
[ -f $O/joint/DONE.json ] || nice -n 10 taskset -c 25-31 $PY m37_joint.py --map $O/align/map.npz --out $O/joint --lambdas 0,0.25,0.5 --workers 4 --cpus 25-31 > $L/m37-joint.out 2> $L/m37-joint.err
echo "joint rc=$? $(date -Is)" >> $O/stage2.log
bash $S/tests/run_tests.sh > $O/tests.log 2>&1
echo "tests rc=$? $(date -Is)" >> $O/stage2.log
nice -n 10 taskset -c 25-31 $PY m37_analysis.py --out $O/results.json > $L/m37-analysis.out 2> $L/m37-analysis.err
echo "analysis rc=$? $(date -Is)" >> $O/stage2.log
touch $O/stage2.done
