#!/bin/bash
# usage: m36_infer.sh <outdir> <model.pt> <regions> [decode_workers=4]
# Retrying, resumable (per region and strand) inference with the P2b wrapper `m34_p2b_infer.py` (m34_infer + M30 registration, so the Dryas / Esch contigs and their
# __val20M pseudo-contigs resolve, also in the spawned decode workers); pc2 backend; the potentials pot/ are kept with M36_KEEP_POT=1 (default; --keep-head).  Waits while the GPU
# has < 5 GB free (the TE project shares it).  Exit 0 when <outdir>/DONE_*.json exists, 9 after 60 failed attempts.
OUT=$1; MODEL=$2; REGIONS=$3; DW=${4:-4}
export M34_BACKEND=pc2
PY=${M36_PY_TRAIN:-$HOME/envs/plantcad2/bin/python}
INFER=${M36_INFER_SCRIPT:-m34_p2b_infer.py}
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
LOGS=${M36_LOGS:-$HOME/coding-rna/logs}
SLEEP_BUSY=${M36_SLEEP_BUSY:-120}; SLEEP_RETRY=${M36_SLEEP_RETRY:-30}
mkdir -p "$OUT" $LOGS
cd ${M36_M34_SCRIPTS:-$HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts}
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
KEEP=""; [ "${M36_KEEP_POT:-1}" = 1 ] && KEEP="--keep-head"
NAME=m36-infer-$(basename $OUT)
free_mib() { if [ -n "$M36_FAKE_FREE_MIB" ]; then echo $M36_FAKE_FREE_MIB; else nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; fi; }
for a in $(seq 1 60); do
  free=$(free_mib)
  if [ "${free:-0}" -lt 5000 ]; then sleep $SLEEP_BUSY; continue; fi
  taskset -c $CPUS $PY $INFER --model $MODEL --regions $REGIONS --out $OUT --decode-workers $DW --cpus $CPUS $KEEP >> $LOGS/$NAME.out 2>> $LOGS/$NAME.err < /dev/null
  rc=$?
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then echo "{\"event\":\"done\",\"time\":\"$(date -Is)\"}" >> $OUT/interruptions.jsonl; exit 0; fi
  echo "{\"event\":\"interruption\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\",\"stderr_tail\":\"$(tail -n 2 $LOGS/$NAME.err | tr '\n"\\' '   ' | cut -c1-300)\"}" >> $OUT/interruptions.jsonl
  sleep $SLEEP_RETRY
done
exit 9
