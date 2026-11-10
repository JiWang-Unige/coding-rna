#!/bin/bash
# PANEL-P6 GPU step: Carbon-A on the whole nuclear genome of every core-10 species, smallest first.  Same pipeline, revision and settings as the M29 G0 (DEV) and M30
# (set B) runs (m30_carbon_raw.sh / m30_native.sh): carbon_a_pipeline_copy model_inference.py with ~/coding-rna/refs/carbon-a-local (HF HuggingFaceBio/Carbon-A-1.2B
# f6789ee), --batch_size 4 --dispatch_chunk_size 1 --attn_implementation sdpa (FlashAttention-2 is not installed; about -0.2 F1 vs the authors' default on one TRAIN
# chromosome, report deck D16), then CPU post_processing.py --viterbi_backend numba --device cpu -> prediction.gff (all genes, confidence= on each gene; the filtered working
# point keeps confidence >= 0.1 or NA at scoring time).
# One input file per main sequence + unplaced.fasta (the other nuclear sequences) per species (p6_genomes.py split; upper case, as in M30).  The post-processing of a file runs
# in the background on cores $PCPUS (at most P6_POST_JOBS at a time) while the GPU moves on; a file's raw_probs.parquet is deleted once its prediction.gff exists.
# Resumable per file; $P6/carbonA/<key>/DONE when every file has its prediction.gff.  usage: p6_carbona.sh [KEY ...]   (P6_SMOKE=1: 2-Mb Arabidopsis training slice)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HERE/p6_lib.sh
KEYS=${*:-$KEYS_BY_SIZE}
NEED=${P6_NEED_MB:-13000}
CA_SRC=$HOME/coding-rna-reorg/experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/src
MODEL=$HOME/coding-rna/refs/carbon-a-local
export PYTHONPATH=$HOME/envs/m29_extra HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
cd $CA_SRC || exit 1
TO=$(( ${P6_TIMEOUT_H:-12} * 3600 ))
POSTJ=${P6_POST_JOBS:-3}

post() {   # KEY STEM OUTDIR : CPU post-processing of one file (background)
  local K=$1 C=$2 O=$3 PQ
  PQ=$(ls $O/*.raw_probs.parquet 2>/dev/null | head -1)
  [ -n "$PQ" ] || { ev event=post_failed key=$K file=$C reason=no_parquet; return 1; }
  mkdir -p $P6/steps
  $GUARD --max-rss-gb ${P6_POST_RSS_GB:-40} --timeout-s 21600 --interval 15 --out $P6/steps/carbon_post_${K}_${C}.res.json -- \
      nice -n 10 taskset -c $PCPUS $PYG post_processing.py --input $PQ --output $O/prediction --viterbi_backend numba --device cpu \
      >> $LOGS/p6-carbon-post-$K.out 2>> $LOGS/p6-carbon-post-$K.err < /dev/null
  local rc=$?
  if [ $rc -eq 0 ] && [ -s $O/prediction.gff ]; then
    echo "{\"parquet\":\"$(basename $PQ)\",\"parquet_bytes\":$(fsize $PQ),\"gff_sha256\":\"$(sha256sum $O/prediction.gff | cut -d' ' -f1)\"}" > $O/post.json
    rm -f $PQ
    $EVT --merge $P6/steps/carbon_post_${K}_${C}.res.json event=post_done key=$K file=$C
  else
    $EVT --merge $P6/steps/carbon_post_${K}_${C}.res.json event=post_failed key=$K file=$C rc=$rc
  fi
}
wait_slots() { while [ "$(jobs -rp | wc -l)" -ge $1 ]; do sleep 20; done; }

infer() {   # KEY STEM FASTA OUTDIR : GPU inference of one file (up to 3 attempts) -> parquet
  local K=$1 C=$2 F=$3 O=$4 a
  for a in 1 2 3 4 5; do                          # 2026-11-10: sporadic SIGSEGV after inference / a tokenizer TypeError, each cured by a plain retry
    gate $NEED carbon_${K}_${C} || return 3
    rm -rf $O; mkdir -p $O
    run_step carbon_${K}_${C} $TO ${P6_MAX_RSS_GB:-48} $PYG model_inference.py --input $F --output_path $O --model_name $MODEL --gpu_count 1 --batch_size 4 \
        --dispatch_chunk_size 1 --attn_implementation sdpa
    ls $O/*.raw_probs.parquet >/dev/null 2>&1 && return 0
    ev event=step_failed step=carbon_${K}_${C} attempt=$a
    sleep 60
  done
  return 1
}

species() {   # KEY DIR INPUTDIR
  local K=$1 D=$2 IN=$3 F C O
  [ -e $D/DONE ] && { ev event=step_skipped step=carbon_$K; return 0; }
  ev event=carbon_species_start key=$K
  for F in $(ls -S $IN/*.fasta); do                # largest file first (keeps the GPU busy while the long post-processings run)
    C=$(basename $F .fasta); O=$D/$C
    [ -s $O/prediction.gff ] && continue
    if ! ls $O/*.raw_probs.parquet >/dev/null 2>&1; then infer $K $C $F $O || return $?; fi
    wait_slots $POSTJ
    post $K $C $O &
  done
  wait
  for F in $IN/*.fasta; do [ -s $D/$(basename $F .fasta)/prediction.gff ] || { ev event=carbon_species_incomplete key=$K file=$(basename $F); return 1; }; done
  touch $D/DONE
  ev event=carbon_species_done key=$K
}

if [ -n "${P6_SMOKE:-}" ]; then
  S=$P6/smoke; mkdir -p $S/carbon_in
  [ -s $S/slice.upper.fa ] || { echo "run P6_SMOKE=1 p6_native.sh first (makes the slice)"; exit 2; }
  cp $S/slice.upper.fa $S/carbon_in/slice.fasta
  species smoke $S/CarbonA $S/carbon_in; exit $?
fi
ev event=tool_start tool=CarbonA keys="$KEYS"
for K in $KEYS; do
  IN=$P6/carbonA/$K/inputs
  if [ ! -e $IN/SPLIT_OK ]; then
    rm -rf $IN; python3.12 -I $P6S/p6_genomes.py split $K $IN > /dev/null && touch $IN/SPLIT_OK || { ev event=tool_aborted tool=CarbonA key=$K reason=split; exit 1; }
  fi
  species $K $P6/carbonA/$K $IN || { rc=$?; ev event=tool_aborted tool=CarbonA key=$K rc=$rc; exit $rc; }
  rm -rf $IN                                        # the inputs are a copy of $G/<key>.nuclear.upper.fa
done
ev event=tool_done tool=CarbonA
exit 0
