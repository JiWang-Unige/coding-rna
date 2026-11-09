#!/bin/bash
# PANEL-P6 GPU step: one native baseline (ANNEVO | Tiberius | Helixer) on the whole nuclear genome of every core-10 species, smallest first.
# Same runtime assets and command lines as the M28 DEV native runs (benchmarks/native_baselines/reports/native_baseline_preparation.md; pipeline/m28/scripts/run_native_baseline.py),
# copied from Baobab outputs/M28-NATIVE-BASELINE-PREP-R2 and M28-NATIVE-RUNTIME208 to $P6/assets:
#   ANNEVO    v2.3.3 @42c920f, ANNEVO_Magnoliopsida.pt, annotation.py -l Magnoliopsida --batch_size 32 -t 8 --num_workers 8 (defaults otherwise; v2.3.3 has no boundary-aware flags)
#             env ~/envs/annevo_p6 = the Baobab `annevo` env pins (python 3.10, torch 2.1.0+cu121, numpy 1.26.4, h5py 3.14.0, biopython 1.87, pandas 2.3.1, numba 0.65.1)
#   Tiberius  2.0.8 source @871d497 (isolated site dir) in the 2.0.7 base SIF, angiosperms.yaml / angiosperms weights, --batch_size 4 --seq_len 200070 --no_softmasking
#   Helixer   0.3.7 SIF (d17bb496), land_plant_v0.3_a_0080.h5, --subsequence-length 64152 --batch-size 16 --deterministic --no-multiprocess
# Differences to DEV (no effect on the predictions): the input is the upper-cased whole nuclear genome ($G/<key>.nuclear.upper.fa), the container bind points are named
# /PANEL-P6-* (so that the GPU gate recognises our container processes), TF_FORCE_GPU_ALLOW_GROWTH=true (TF does not grab the whole card), and on an out-of-memory failure
# the batch size is halved (recorded in run.json).
# Completion = the native completion marker of the M28 runner in this step's logs + p6_check_gff.py; then $P6/native/<tool>/<key>/DONE.  Resumable per genome.
# usage: p6_native.sh TOOL [KEY ...]        (P6_SMOKE=1: a 2-Mb Arabidopsis training slice instead, not scored)
TOOL=$1; shift
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HERE/p6_lib.sh
PR=$A/M28-NATIVE-BASELINE-PREP-R2
RT=$A/M28-NATIVE-RUNTIME208
KEYS=${*:-$KEYS_BY_SIZE}
case $TOOL in
  ANNEVO)   B0=32; NEED=${P6_NEED_MB:-10000}; MARK='The gene annotation took ';;
  Tiberius) B0=4;  NEED=${P6_NEED_MB:-12000}; MARK='Tiberius took ';;
  Helixer)  B0=16; NEED=${P6_NEED_MB:-12000}; MARK='Helixer successfully finished the annotation';;
  *) echo "unknown tool $TOOL"; exit 2;;
esac
TO=$(( ${P6_TIMEOUT_H:-30} * 3600 ))
TFENV=(--env CUDA_VISIBLE_DEVICES=0 --env OMP_NUM_THREADS=8 --env TF_NUM_INTRAOP_THREADS=8 --env TF_NUM_INTEROP_THREADS=1 --env TF_FORCE_GPU_ALLOW_GROWTH=true)

cmd_for() {   # KEY FASTA_DIR FASTA_NAME RUN_DIR BATCH
  local K=$1 FD=$2 FN=$3 D=$4 B=$5
  case $TOOL in
    ANNEVO) echo $HOME/envs/annevo_p6/bin/python $PR/assets/ANNEVO/annotation.py -g $FD/$FN -m $PR/assets/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt -l Magnoliopsida \
              -o $D/prediction.gff3 --batch_size $B -t 8 --num_workers 8 --tmp_path $D/temporary --show_log;;
    Tiberius) echo $APPT exec --nv --cleanenv "${TFENV[@]}" --bind $PR:/m28-prep:ro --bind $D:/PANEL-P6-run --pwd /PANEL-P6-run --bind $RT:/m28-runtime:ro \
              --env PYTHONPATH=/m28-runtime/site --bind $PR/assets/Tiberius-model-weights:/m28-runtime/site/model_weights:ro --bind $FD:/PANEL-P6-genomes:ro \
              $PR/assets/Tiberius-2.0.8.sif python3 -m tiberius.main --genome /PANEL-P6-genomes/$FN --model_cfg /m28-prep/assets/angiosperms.yaml \
              --out /PANEL-P6-run/prediction.gff3 --batch_size $B --seq_len 200070 --no_softmasking;;
    Helixer) echo $APPT exec --nv --cleanenv "${TFENV[@]}" --bind $PR:/m28-prep:ro --bind $D:/PANEL-P6-run --pwd /home/helixer_user/Helixer --bind $FD:/PANEL-P6-genomes:ro \
              $PR/assets/Helixer-0.3.7.sif python3 /home/helixer_user/Helixer/Helixer.py --fasta-path /PANEL-P6-genomes/$FN --gff-output-path /PANEL-P6-run/prediction.gff3 \
              --species $K --model-filepath /m28-prep/assets/land_plant_v0.3_a_0080.h5 --subsequence-length 64152 --batch-size $B --temporary-dir /PANEL-P6-run/temporary \
              --deterministic --no-multiprocess;;
  esac
}
oom_since() { tail -c +$(( $2 + 1 )) "$1" 2>/dev/null | grep -qiE 'ResourceExhausted|out of memory|OOM when allocating|CUDA_ERROR_OUT_OF_MEMORY|CUDA out of memory'; }

one() {   # KEY FASTA_DIR FASTA_NAME RUN_DIR
  local K=$1 FD=$2 FN=$3 D=$4 NAME=${TOOL}_$1 B=$B0 rc o0 e0 attempt chk
  [ -e $D/DONE ] && { ev event=step_skipped step=$NAME; return 0; }
  for attempt in 1 2 3 4; do
    gate $NEED $NAME || return 3
    rm -rf $D; mkdir -p $D/temporary                    # ANNEVO needs --tmp_path to exist (the M28 DEV runner pre-created it)
    o0=$(fsize $LOGS/p6-$NAME.out); e0=$(fsize $LOGS/p6-$NAME.err)
    run_step $NAME $TO ${P6_MAX_RSS_GB:-64} $(cmd_for $K $FD $FN $D $B)
    rc=$?
    if [ $rc -eq 0 ] && { tail -c +$((o0+1)) $LOGS/p6-$NAME.out; tail -c +$((e0+1)) $LOGS/p6-$NAME.err; } | grep -qF "$MARK" \
       && chk=$(python3.12 -I $HERE/p6_check_gff.py $D/prediction.gff3 $FD/$FN ${LENJSON:+$LENJSON} 2>&1); then
      python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_event.py --file $D/run.json --write --merge $P6/steps/$NAME.res.json tool=$TOOL key=$K \
          batch=$B "fasta=$FD/$FN" "check=$chk" "gff_sha256=$(sha256sum $D/prediction.gff3 | cut -d' ' -f1)" "gpu=$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader)"
      rm -rf $D/temporary
      touch $D/DONE
      ev event=step_done step=$NAME batch=$B "check=$chk"
      return 0
    fi
    ev event=step_failed step=$NAME attempt=$attempt rc=$rc "check=${chk:-}"
    if oom_since $LOGS/p6-$NAME.err $e0 || oom_since $LOGS/p6-$NAME.out $o0; then
      [ $B -le 1 ] && return 1
      B=$(( B / 2 )); ev event=oom_retry step=$NAME batch=$B
    else
      [ $rc -eq 137 ] && return 1                     # killed by the RAM / time guard: never retried
      [ $attempt -ge 2 ] && return 1
      sleep 60
    fi
  done
  return 1
}

if [ -n "${P6_SMOKE:-}" ]; then                     # end-to-end smoke on a 2-Mb Arabidopsis TRAINING slice (M38 smoke region), never scored
  S=$P6/smoke; mkdir -p $S
  [ -s $S/slice.upper.fa ] || { python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_slice.py $HOME/coding-rna/outputs/M38-VIPSANIA/genomes/arabidopsis_thaliana.nuclear.fa NC_003070.9 5000001 7000000 $S/slice.fa >/dev/null && awk '/^>/{print;next}{print toupper($0)}' $S/slice.fa > $S/slice.upper.fa; }
  NEED=${P6_NEED_MB:-$NEED} TO=1800 one smoke $S slice.upper.fa $S/$TOOL
  exit $?
fi
ev event=tool_start tool=$TOOL keys="$KEYS"
for K in $KEYS; do
  LENJSON=$HOME/coding-rna/data/panel/prepared/$K/lengths.json one $K $G $K.nuclear.upper.fa $P6/native/$TOOL/$K || { rc=$?; ev event=tool_aborted tool=$TOOL key=$K rc=$rc; exit $rc; }
done
ev event=tool_done tool=$TOOL
exit 0
