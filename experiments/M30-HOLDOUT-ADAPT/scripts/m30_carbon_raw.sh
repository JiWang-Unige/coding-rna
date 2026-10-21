#!/bin/bash
# Carbon-A native raw-probability inference (SDPA, as in M29 G0) for every M30 contig; resumable per contig; one GPU job at a time.
O=$HOME/coding-rna/outputs/M30-HOLDOUT-ADAPT/carbonA_raw
cd $HOME/coding-rna-reorg/experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/src
export PYTHONPATH=$HOME/envs/m29_extra HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
for SP in dryas_octopetala eschscholzia_californica; do
  for F in $HOME/coding-rna/data/m1_screen/$SP/inputs/*.fasta; do
    C=$(basename $F .fasta); D=$O/$SP/$C
    if ls $D/*.raw_probs.parquet >/dev/null 2>&1; then continue; fi
    for a in 1 2 3; do
      free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
      echo "$(date -Is) start $SP $C attempt $a free=$free" >> $O/log.txt
      rm -rf $D; mkdir -p $D
      taskset -c 2,4,6,8,10,12,14,16-24 $HOME/envs/generanno/bin/python model_inference.py --input $F --output_path $D --model_name $HOME/coding-rna/refs/carbon-a-local --gpu_count 1 --batch_size 4 --dispatch_chunk_size 1 --attn_implementation sdpa > $HOME/coding-rna/logs/m30-carbonraw-$C.out 2> $HOME/coding-rna/logs/m30-carbonraw-$C.err < /dev/null
      if ls $D/*.raw_probs.parquet >/dev/null 2>&1; then echo "$(date -Is) done $C" >> $O/log.txt; break; fi
      echo "$(date -Is) FAILED $C attempt $a: $(tail -n 2 $HOME/coding-rna/logs/m30-carbonraw-$C.err | tr '\n' ' ' | cut -c1-200)" >> $O/log.txt; sleep 20
    done
  done
done
echo "$(date -Is) ALL_DONE" >> $O/log.txt
