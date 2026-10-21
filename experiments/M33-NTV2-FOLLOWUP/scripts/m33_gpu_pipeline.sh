#!/bin/bash
# GPU stage of M33 (run ONLY after M32_DONE exists and no M32 process is left; the caller checks).  Sequential, resumable per (region, strand):
#   val (internal-validation chromosomes) -> dev -> Dryas B -> Eschscholzia B.     One GPU job at a time; retries when the GPU is full (other users).
# usage: m33_gpu_pipeline.sh            logs: ~/coding-rna/logs/m33-infer-<name>.{out,err}; marker: $O/gpu_inference.done
O=$HOME/coding-rna/outputs/M33-NTV2-FOLLOWUP
MODEL=$O/model_B_step1536.pt
CPUS=2,4,6,8,10,12,14,16
cd $HOME/coding-rna-reorg/experiments/M33-NTV2-FOLLOWUP/scripts
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
run() {   # name regions
  NAME=$1; REG=$2; OUT=$O/$NAME; mkdir -p $OUT
  for a in $(seq 1 30); do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "$free" -lt 5000 ]; then sleep 120; continue; fi
    taskset -c $CPUS $HOME/envs/generanno/bin/python m33_infer.py --model $MODEL --regions $REG --out $OUT --decode-workers 4 --cpus $CPUS >> $HOME/coding-rna/logs/m33-infer-$NAME.out 2>> $HOME/coding-rna/logs/m33-infer-$NAME.err < /dev/null
    rc=$?
    if ls $OUT/DONE_*.json >/dev/null 2>&1; then echo "{\"event\":\"done\",\"name\":\"$NAME\",\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl; return 0; fi
    echo "{\"event\":\"interruption\",\"name\":\"$NAME\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\",\"stderr_tail\":\"$(tail -n 2 $HOME/coding-rna/logs/m33-infer-$NAME.err | tr '\n"' '  ' | cut -c1-300)\"}" >> $O/interruptions.jsonl
    sleep 30
  done
  return 1
}
run val NC_003071.7,NC_089044.1 || exit 1
run dev NC_003074.8,NC_089041.1 || exit 1
run evalB_dryas NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1 || exit 1
run evalB_esch NW_028238309.1,NW_028238461.1,NW_028238483.1 || exit 1
echo "{\"event\":\"gpu_inference_done\",\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl
touch $O/gpu_inference.done
