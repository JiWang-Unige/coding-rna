#!/bin/bash
# PANEL-P6 GPU step: Vipsania on the whole nuclear genome (NCBI soft-masking kept) of every core-10 species, smallest first.  Same environment, model and settings as M38
# (report 30 §1-3; experiments/M38-VIPSANIA/scripts/vip_devb.sh): ~/envs/vipsania_m38 (vipsania 1.1.0), Streptophyta model pinned by ID j9m0cdmk (sha256 checked),
#   ft  (C1-criterion mode, report 28 §6): vipsania annotate j9m0cdmk <genome> --finetune  (author defaults: 10 epochs, lr 1e-4, effective finetuning batch 64)
#   pre (descriptive):                     vipsania annotate j9m0cdmk <genome>
# Only the settings forced by the 24-GB card differ (as in M38): --finetune_B 2 when >= 23,040 MiB are free (the step first waits up to P6_VIP_WAIT_B2_H = 2 h for that),
# else 1 (the effective batch stays 64 by gradient accumulation); after an out-of-memory failure a smaller annotation batch -B.  No parameter is tuned, nothing is re-run
# because of a result; only failures are retried.
# Dryas / Eschscholzia pre: the M38 set-B run already annotated the identical genome (sha256 equal, genomes_manifest.json same_as_m38_nuclear_fa) with the same model and env;
# its GFF is copied instead of re-running (recorded in run.json as reused_from).
# Resumable per genome and mode: $P6/vipsania/runs/<key>/<mode>/<mode>_DONE.   usage: p6_vipsania.sh [KEY ...]   (P6_SMOKE=1: pre on a 2-Mb Arabidopsis training slice)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HERE/p6_lib.sh
M38=$HOME/coding-rna/outputs/M38-VIPSANIA
VENV=$HOME/envs/vipsania_m38
VIP=$VENV/bin/vipsania
export VIPSANIA_CACHE=$M38/vipsania_cache MODEL_ID=j9m0cdmk
export TF_CPP_MIN_LOG_LEVEL=2 WANDB_MODE=offline WANDB_SILENT=true WANDB_DISABLED=true TF_FORCE_GPU_ALLOW_GROWTH=true
KEYS=${*:-$KEYS_BY_SIZE}
FT_TO=$(( ${P6_VIP_FT_TIMEOUT_H:-16} * 3600 )); ANN_TO=$(( ${P6_VIP_ANN_TIMEOUT_H:-5} * 3600 ))
MIN_GENES=${P6_VIP_MIN_GENES:-1000}
(cd $VIPSANIA_CACHE/models/$MODEL_ID && sha256sum -c --quiet $M38/model_$MODEL_ID.sha256) >/dev/null 2>&1 || { ev event=refused reason="model $MODEL_ID missing or checksum differs"; exit 11; }
export LD_LIBRARY_PATH=$($VENV/bin/python -c "import site, glob, os; print(':'.join(sorted(glob.glob(os.path.join(site.getsitepackages()[0], 'nvidia', '*', 'lib')))))" 2>/dev/null)${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
oom_since() { tail -c +$(( $2 + 1 )) "$1" 2>/dev/null | grep -qiE 'ResourceExhausted|out of memory|OOM when allocating|CUDA_ERROR_OUT_OF_MEMORY|Allocator .* ran out'; }
genes_in() { awk -F'\t' '$3=="gene"{n++} END{print n+0}' "$1" 2>/dev/null; }
log_field() { grep -m1 "$2" "$1" 2>/dev/null | sed 's/^[> ]*//' | tr -d '"\\' ; }

wait_b2() {   # wait up to P6_VIP_WAIT_B2_H for >= 23,040 MiB free (finetune_B 2); 0 if reached
  local left=$(( ${P6_VIP_WAIT_B2_H:-2} * 3600 )) f
  while [ $left -gt 0 ]; do f=$(gpu_free); [ "${f:-0}" -ge 23040 ] && [ -z "$(gpu_blockers)" ] && return 0; sleep 120; left=$((left-120)); done
  return 1
}

annotate() {   # KEY MODE FASTA RUNDIR
  local K=$1 MODE=$2 FA=$3 D=$4 NAME=vip_${1}_$2 TO rc e0 B=0 ann_b=() extra=() attempt f
  [ -e $D/${MODE}_DONE ] && { ev event=step_skipped step=$NAME; return 0; }
  mkdir -p $D
  TO=$ANN_TO; [ $MODE = ft ] && TO=$FT_TO
  for attempt in 1 2 3 4; do
    gate 11520 $NAME || return 3
    extra=()
    if [ $MODE = ft ]; then
      [ "$B" = 1 ] || { wait_b2; gate 11520 $NAME || return 3; }
      f=$(gpu_free); B=1; [ "${f:-0}" -ge 23040 ] && [ "${FORCE_B1:-0}" = 0 ] && B=2
      extra=(--finetune --finetune_B $B)
    fi
    rm -rf $D/vipsania.gff3 $D/vipsania.log $D/finetuning_*
    e0=$(fsize $LOGS/p6-$NAME.err)
    run_step $NAME $TO ${P6_VIP_MAX_RSS_GB:-48} $VIP annotate $MODEL_ID $FA -o $D/vipsania.gff3 "${ann_b[@]}" "${extra[@]}"
    rc=$?
    if [ $rc -eq 0 ] && [ "$(genes_in $D/vipsania.gff3)" -ge $MIN_GENES ]; then
      python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_event.py --file $D/run.json --write --merge $P6/steps/$NAME.res.json genome=$K mode=$MODE model=$MODEL_ID \
          "fasta=$FA" genes=$(genes_in $D/vipsania.gff3) "finetune_B=$B" "annotate_B_flag=${ann_b[*]:-auto}" "batch_log=$(log_field $D/vipsania.log 'batch size for maximal chunk length')" \
          "peak_gpu_log=$(log_field $D/vipsania.log 'peak GPU memory')" "log_time=$(log_field $D/vipsania.log '> time:')" "gff_sha256=$(sha256sum $D/vipsania.gff3 | cut -d' ' -f1)"
      touch $D/${MODE}_DONE
      ev event=step_done step=$NAME genes=$(genes_in $D/vipsania.gff3) finetune_B=$B
      return 0
    fi
    ev event=step_failed step=$NAME attempt=$attempt rc=$rc genes=$(genes_in $D/vipsania.gff3)
    if oom_since $LOGS/p6-$NAME.err $e0; then
      if [ $MODE = ft ] && [ "$B" -gt 1 ]; then FORCE_B1=1; B=1; ev event=oom_retry step=$NAME finetune_B=1
      elif [ ${#ann_b[@]} -eq 0 ]; then ann_b=(-B 4); ev event=oom_retry step=$NAME annotate_B=4
      elif [ "${ann_b[1]}" -gt 1 ]; then ann_b=(-B $(( ${ann_b[1]} / 2 ))); ev event=oom_retry step=$NAME annotate_B=${ann_b[1]}
      else return 1; fi
    else
      [ $rc -eq 137 ] && return 1
      [ $attempt -ge 2 ] && return 1
      sleep 60
    fi
  done
  return 1
}

reuse_m38_pre() {   # KEY RUNDIR : copy the M38 pre annotation of the identical genome
  local K=$1 D=$2
  python3.12 -I - <<PY || return 1
import json, sys
m = json.load(open("$G/genomes_manifest.json"))["$K"]
sys.exit(0 if m.get("same_as_m38_nuclear_fa") else 1)
PY
  [ -e $M38/runs/$K/pre/pre_DONE ] || return 1
  mkdir -p $D && cp $M38/runs/$K/pre/vipsania.gff3 $D/vipsania.gff3 && cp $M38/runs/$K/pre/run.json $D/run_m38.json || return 1
  python3.12 -I $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_event.py --file $D/run.json --write genome=$K mode=pre model=$MODEL_ID \
      "reused_from=$M38/runs/$K/pre" genes=$(genes_in $D/vipsania.gff3) "gff_sha256=$(sha256sum $D/vipsania.gff3 | cut -d' ' -f1)"
  touch $D/pre_DONE
  ev event=step_reused step=vip_${K}_pre from=$M38/runs/$K/pre
}

if [ -n "${P6_SMOKE:-}" ]; then
  S=$P6/smoke; [ -s $S/slice.fa ] || { echo "run P6_SMOKE=1 p6_native.sh first (makes the slice)"; exit 2; }
  MIN_GENES=100 annotate smoke pre $S/slice.fa $S/Vipsania_pre; exit $?
fi
ev event=tool_start tool=Vipsania keys="$KEYS" model=$MODEL_ID
for K in $KEYS; do
  R=$P6/vipsania/runs/$K
  annotate $K ft $G/$K.nuclear.fa $R/ft || { rc=$?; ev event=tool_aborted tool=Vipsania step=${K}_ft rc=$rc; exit $rc; }
  [ -e $R/pre/pre_DONE ] || reuse_m38_pre $K $R/pre || annotate $K pre $G/$K.nuclear.fa $R/pre || { rc=$?; ev event=tool_aborted tool=Vipsania step=${K}_pre rc=$rc; exit $rc; }
done
ev event=tool_done tool=Vipsania
exit 0
