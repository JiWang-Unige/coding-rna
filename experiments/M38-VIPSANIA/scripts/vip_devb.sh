#!/bin/bash
# M38-VIPSANIA GPU step (queue m34_after_t2s2.sh, step "vipsania"; user decision vipsania-run 2026-11-04).  Protocol: docs/review_20261016/28_animal_C2_and_Vipsania_proposal.md section 6.
#   0 gate      refuses to start without $M/READY (written by hand after review), exits quietly if $M/SKIP or $M/DEVB_DONE exists; single instance (flock)
#   1 GPU smoke pretrained annotate of a 2-Mb TRAINING slice of Arabidopsis NC_003070.9 on the GPU (must report a GPU, else exit 20) + a short finetune probe that measures the
#               finetune batch size that fits and the seconds per window (-> projected finetune hours; exit 25 if the upper-bound projection exceeds VIP_MAX_PROJ_H = 6 h per genome; user approved ~15-20 GPU h for DEV + B, 2026-11-04).  Not scored.
#   2 genomes   smallest first: arabidopsis_thaliana (DEV) 119 Mb, dryas_octopetala (B) 230 Mb, oryza_sativa (DEV) 386 Mb, eschscholzia_californica (B) 401 Mb; for each the
#               whole nuclear genome is annotated twice: pretrained only (pre) and `--finetune` (ft, the C1-criterion mode), author defaults (10 epochs, lr 1e-4, T 200,000,
#               effective finetuning batch 64), model pinned by ID j9m0cdmk.  Only settings forced by the 24-GB card differ: --finetune_B (2 if >= 23,040 MiB are free, else 1; the
#               effective batch stays 64 through gradient accumulation) and, after an out-of-memory failure, a smaller annotation batch -B.
#   3 scoring   vip_score.py dev / setb (CPU) on both modes, then $M/DEVB_DONE.
# Resumable: step sentinels $M/runs/<genome>/{pre,ft}/*_DONE and the per-genome sentinel $M/runs/<genome>/DONE are skipped on a restart (a half-finished finetuning restarts).
# Events (one JSON line per start/end): $M/events.jsonl ; logs: ~/coding-rna/logs/m38-vipsania-<step>.{out,err} ; per-step resource files $M/steps/<step>.res.json.
# Status:  tail -n 5 ~/coding-rna/outputs/M38-VIPSANIA/events.jsonl
# Test (no GPU, no Vipsania): tests/test_devb_stub.sh  (stub binary, throw-away M38_ROOT; the real READY file is never touched)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/vip_env.sh"
M=$M38
LOGS=${VIP_LOGS:-$HOME/coding-rna/logs}
EV=$M/events.jsonl
CPUS=${VIP_CPUS:-2,4,6,8,10,12,14,16-24}               # the cores the M34 GPU steps use (the TE project's CPU jobs keep the others)
SCPUS=${VIP_SCORE_CPUS:-18-23}
FT_TIMEOUT_S=$(( ${VIP_FT_TIMEOUT_H:-16} * 3600 ))
ANN_TIMEOUT_S=$(( ${VIP_ANN_TIMEOUT_H:-5} * 3600 ))
MIN_GENES=${VIP_MIN_GENES:-1000}
GENOMES="arabidopsis_thaliana dryas_octopetala oryza_sativa eschscholzia_californica"      # ascending genome size
DEV_GENOMES="arabidopsis_thaliana oryza_sativa"
B_GENOMES="dryas_octopetala eschscholzia_californica"
FT_GENOMES=${VIP_FT_GENOMES:-$GENOMES}   # genomes that also get the --finetune run (user decision vipsania-ft-budget 2026-11-04: DEV only; set B finetune comes from the panel run)
in_ft() { case " $FT_GENOMES " in *" $1 "*) return 0;; esac; return 1; }
GUARD="python3.12 -I $HERE/vip_guard.py"
EVT="python3.12 -I $HERE/vip_event.py --file $EV"
mkdir -p $M/steps $M/runs $LOGS
ev() { $EVT "$@"; }

# ---- 0 gate ---------------------------------------------------------------------------------------------------------------------------------------------------
[ -e $M/SKIP ] && { ev event=skipped reason=SKIP; exit 0; }
[ -e $M/DEVB_DONE ] && { ev event=already_done; exit 0; }
[ -e $M/READY ] || { ev event=refused reason=no_READY; exit 11; }
exec 9>$M/.devb.lock
flock -n 9 || { ev event=refused reason=already_running; exit 3; }
[ -x $VIP ] || { ev event=refused reason="no $VIP"; exit 11; }
[ "${CUDA_VISIBLE_DEVICES+x}" = x ] && [ -z "$CUDA_VISIBLE_DEVICES" ] && unset CUDA_VISIBLE_DEVICES      # an empty value would hide the GPU
(cd $VIPSANIA_CACHE/models/$MODEL_ID && sha256sum -c --quiet $M/model_$MODEL_ID.sha256) >/dev/null 2>&1 || { ev event=refused reason="model $MODEL_ID missing or checksum differs"; exit 11; }
for g in $GENOMES; do
  [ -s $M/genomes/$g.nuclear.fa ] || { ev event=refused reason="genome $g missing"; exit 11; }
done
python3.12 -I - <<PY || { ev event=refused reason="genome FASTA size differs from genomes_manifest.json"; exit 11; }
import json, os, sys
m = json.load(open("$M/genomes/genomes_manifest.json"))
sys.exit(1 if [k for k in "$GENOMES".split() if os.path.getsize(f"$M/genomes/{k}.nuclear.fa") != m[k]["bytes"]] else 0)
PY
export TF_FORCE_GPU_ALLOW_GROWTH=true
export LD_LIBRARY_PATH=$($VPY -c "import site, glob, os; print(':'.join(sorted(glob.glob(os.path.join(site.getsitepackages()[0], 'nvidia', '*', 'lib')))))" 2>/dev/null)${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
ev event=devb_start model=$MODEL_ID genomes="$GENOMES" cpus=$CPUS "gpu=$(nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader | tr -d ' ')"

gpu_free() { [ -n "${VIP_GPU_FREE_MB:-}" ] && { echo $VIP_GPU_FREE_MB; return; }   # (override for the stub test only)
  nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -dc 0-9; }
wait_free() {   # min_MiB max_wait_s : wait until the GPU has that much free memory (2 consecutive checks 10 s apart); 1 on timeout
  local need=$1 left=$2 ok=0 f
  while :; do
    f=$(gpu_free); f=${f:-0}
    if [ "$f" -ge "$need" ]; then ok=$((ok+1)); [ $ok -ge 2 ] && return 0; sleep 10
    else ok=0; left=$((left-60)); [ $left -le 0 ] && return 1; sleep 60; fi
  done
}
pick_ftb() { local f; f=$(gpu_free); f=${f:-0}; if [ "$f" -ge 23040 ]; then echo 2; elif [ "$f" -ge 11520 ]; then echo 1; else echo 0; fi; }   # the authors' own cost model: ~10.1 GiB per 20-kb window, 0.9 safety
oom_since() { tail -c +$(( $2 + 1 )) "$1" 2>/dev/null | grep -qiE 'ResourceExhausted|out of memory|OOM when allocating|CUDA_ERROR_OUT_OF_MEMORY|Allocator .* ran out'; }
fsize() { stat -c %s "$1" 2>/dev/null || echo 0; }

run_step() {   # name timeout_s cmd... : one command under the RAM / timeout / GPU-memory accounting guard; stdout/stderr appended to ~/coding-rna/logs/m38-vipsania-<name>.{out,err}
  local NAME=$1 TO=$2; shift 2
  local OUT=$LOGS/m38-vipsania-$NAME.out ERR=$LOGS/m38-vipsania-$NAME.err RES=$M/steps/$NAME.res.json
  rm -f $RES
  ev event=step_start step=$NAME "gpu_free_mb=$(gpu_free)" "cmd=$*"
  $GUARD --gpu --max-rss-gb ${VIP_MAX_RSS_GB:-48} --timeout-s $TO --interval 10 --out $RES -- nice -n 10 taskset -c $CPUS "$@" >> $OUT 2>> $ERR < /dev/null
  local rc=$?
  $EVT --merge $RES event=step_end step=$NAME
  return $rc
}
genes_in() { awk -F'\t' '$3=="gene"{n++} END{print n+0}' "$1" 2>/dev/null; }
log_field() { grep -m1 "$2" "$1" 2>/dev/null | sed 's/^[> ]*//' | tr -d '"\\' ; }

# ---- 1 GPU smoke ----------------------------------------------------------------------------------------------------------------------------------------------
S=$M/gpu_smoke; mkdir -p $S/pre $S/ft
smoke_gpu() {
  [ -e $S/GPU_SMOKE_OK ] && { ev event=step_skipped step=gpu_smoke; return 0; }
  [ -s $S/slice.fa ] || python3.12 -I $HERE/vip_slice.py $M/genomes/arabidopsis_thaliana.nuclear.fa NC_003070.9 5000001 7000000 $S/slice.fa >/dev/null || return 20
  wait_free 11520 7200 || { ev event=gpu_smoke_failed reason="less than 11.5 GiB of GPU memory free for 2 h"; return 20; }
  local o0 e0 rc
  o0=$(fsize $LOGS/m38-vipsania-gpusmoke_pre.out); e0=$(fsize $LOGS/m38-vipsania-gpusmoke_pre.err)
  run_step gpusmoke_pre $ANN_TIMEOUT_S $VIP annotate $MODEL_ID $S/slice.fa -o $S/pre/smoke.gff3
  rc=$?
  if [ $rc -ne 0 ] || ! tail -c +$(( o0 + 1 )) $LOGS/m38-vipsania-gpusmoke_pre.out | grep -qE '^Using [0-9]+ GPU'; then
    ev event=gpu_smoke_failed reason="pretrained smoke rc=$rc or no GPU reported (see $LOGS/m38-vipsania-gpusmoke_pre.{out,err})"; return 20
  fi
  [ "$(genes_in $S/pre/smoke.gff3)" -ge ${VIP_SMOKE_MIN_GENES:-100} ] || { ev event=gpu_smoke_failed reason="too few genes in the smoke annotation"; return 20; }
  # finetune probe: 512 windows (B x train_steps), effective batch 4 x B only for the probe (author default 64 is used in production)
  local B windows=512 ts
  B=$(pick_ftb)
  while [ "$B" -ge 1 ]; do
    ts=$(( windows / B ))
    o0=$(fsize $LOGS/m38-vipsania-gpusmoke_ft.err); e0=$o0
    run_step gpusmoke_ft $FT_TIMEOUT_S $VIP annotate $MODEL_ID $S/slice.fa -o $S/ft/smoke.gff3 --finetune --finetune_B $B --finetune_epochs 1 \
        -fo trainer.train_steps=$ts -fo trainer.gradient_accumulation_steps=4
    rc=$?
    if [ $rc -eq 0 ]; then break; fi
    if oom_since $LOGS/m38-vipsania-gpusmoke_ft.err $o0 && [ "$B" -gt 1 ]; then ev event=ft_probe_oom finetune_B=$B; B=$((B/2)); else B=0; fi
  done
  [ "$B" -ge 1 ] || { ev event=gpu_smoke_failed reason="finetune probe failed at every batch size"; return 20; }
  echo $B > $M/ft_B.txt
  python3.12 -I - <<PY >> $EV
import json, time
pre = json.load(open("$M/steps/gpusmoke_pre.res.json")); ft = json.load(open("$M/steps/gpusmoke_ft.res.json"))
win = $windows
t_win = max(ft["wall_seconds"] - pre["wall_seconds"], 1) / win           # upper bound: includes the training-step compile
ev = dict(event="gpu_smoke_done", finetune_B=$B, probe_windows=win, pre_seconds=pre["wall_seconds"], ft_probe_seconds=ft["wall_seconds"],
          sec_per_window_upper=round(t_win, 3), projected_finetune_hours_per_genome_upper=round(t_win * 64000 / 3600, 2),
          peak_gpu_mb_pre=pre["peak_gpu_mb"], peak_gpu_mb_ft=ft["peak_gpu_mb"], time=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
open("$S/projection_h.txt", "w").write(str(ev["projected_finetune_hours_per_genome_upper"]))
print(json.dumps(ev))
PY
  touch $S/GPU_SMOKE_OK
  return 0
}

# ---- 2 one genome, one mode -------------------------------------------------------------------------------------------------------------------------------------
annotate() {   # genome mode(pre|ft)
  local G=$1 MODE=$2 D=$M/runs/$1/$2 FA=$M/genomes/$1.nuclear.fa NAME=${1}_$2 TO rc o0 B=0 ann_b=() extra=() attempt
  [ -e $D/${MODE}_DONE ] && { ev event=step_skipped step=$NAME; return 0; }
  mkdir -p $D
  TO=$ANN_TIMEOUT_S; [ $MODE = ft ] && TO=$FT_TIMEOUT_S
  for attempt in 1 2 3 4; do
    wait_free 11520 7200 || { ev event=step_failed step=$NAME reason="no free GPU memory for 2 h"; return 1; }
    extra=(); [ $MODE = ft ] && { B=$(cat $M/ft_B.txt 2>/dev/null || echo 2); f=$(gpu_free); [ "${f:-0}" -lt 23040 ] && B=1; extra=(--finetune --finetune_B $B); }
    rm -rf $D/vipsania.gff3 $D/vipsania.log $D/finetuning_*          # leftovers of a failed attempt (the finetuned checkpoint is ~400 MB)
    o0=$(fsize $LOGS/m38-vipsania-$NAME.err)
    run_step $NAME $TO $VIP annotate $MODEL_ID $FA -o $D/vipsania.gff3 "${ann_b[@]}" "${extra[@]}"
    rc=$?
    if [ $rc -eq 0 ] && [ "$(genes_in $D/vipsania.gff3)" -ge $MIN_GENES ]; then
      python3.12 -I $HERE/vip_event.py --file $D/run.json --write --merge $M/steps/$NAME.res.json genome=$G mode=$MODE model=$MODEL_ID \
          "fasta=$FA" genes=$(genes_in $D/vipsania.gff3) "finetune_B=$B" "annotate_B_flag=${ann_b[*]:-auto}" "batch_log=$(log_field $D/vipsania.log 'batch size for maximal chunk length')" \
          "peak_gpu_log=$(log_field $D/vipsania.log 'peak GPU memory')" "log_time=$(log_field $D/vipsania.log '> time:')" "gff_sha256=$(sha256sum $D/vipsania.gff3 | cut -d' ' -f1)"
      touch $D/${MODE}_DONE
      ev event=step_done step=$NAME genes=$(genes_in $D/vipsania.gff3)
      return 0
    fi
    ev event=step_failed step=$NAME attempt=$attempt rc=$rc
    if oom_since $LOGS/m38-vipsania-$NAME.err $o0; then       # forced by the 24-GB card: smaller batches
      if [ $MODE = ft ] && [ "$B" -gt 1 ]; then echo 1 > $M/ft_B.txt; ev event=oom_retry step=$NAME finetune_B=1
      elif [ ${#ann_b[@]} -eq 0 ]; then ann_b=(-B 4); ev event=oom_retry step=$NAME annotate_B=4
      elif [ "${ann_b[1]}" -gt 1 ]; then ann_b=(-B $(( ${ann_b[1]} / 2 ))); ev event=oom_retry step=$NAME annotate_B=${ann_b[1]}
      else return 1; fi
    else
      [ $rc -eq 137 ] && return 1                      # killed by the RAM / time guard: never retried
      sleep 60
    fi
  done
  return 1
}

# ---- run --------------------------------------------------------------------------------------------------------------------------------------------------------
smoke_gpu || { rc=$?; ev event=devb_aborted stage=gpu_smoke rc=$rc; exit ${rc:-20}; }
# protect the shared queue from a runaway plan: the finetuning is a fixed 64,000 windows per genome; if the measured (upper-bound) projection exceeds VIP_MAX_PROJ_H hours per
# genome the job stops here (exit 25) and the queue moves on to T2_s3 -- the decision to run it anyway (VIP_MAX_PROJ_H=99 or a longer slot) is the user's
PROJ=$(cat $S/projection_h.txt 2>/dev/null || echo 0)
if awk "BEGIN{exit !($PROJ > ${VIP_MAX_PROJ_H:-6})}"; then ev event=devb_aborted stage=projection projected_finetune_hours_per_genome_upper=$PROJ max=${VIP_MAX_PROJ_H:-6}; exit 25; fi
for g in $GENOMES; do
  if [ -e $M/runs/$g/DONE ]; then ev event=genome_skipped genome=$g; continue; fi
  annotate $g pre || { ev event=devb_aborted stage="${g}_pre"; exit 21; }
  if in_ft $g; then annotate $g ft || { ev event=devb_aborted stage="${g}_ft"; exit 22; }; else ev event=step_skipped step=${g}_ft reason=not_in_FT_GENOMES; fi
  echo "{\"genome\":\"$g\",\"time\":\"$(date -Is)\"}" > $M/runs/$g/DONE
  ev event=genome_done genome=$g
done

# ---- 3 scoring (CPU, generanno env: the project's ruler) ----------------------------------------------------------------------------------------------------------
if [ -z "${VIP_SKIP_SCORE:-}" ]; then
  R=$M/runs
  score() {   # kind(dev|setb) genomes...
    local KIND=$1; shift; local ft=() pre=() g
    for g in "$@"; do in_ft $g && ft+=("$g=$R/$g/ft/vipsania.gff3"); pre+=("$g=$R/$g/pre/vipsania.gff3"); done
    local ftarg=(); [ ${#ft[@]} -gt 0 ] && ftarg=(--ft "${ft[@]}")
    [ -e $M/results/$KIND/vip_${KIND}_result.json ] && { ev event=step_skipped step=score_$KIND; return 0; }
    ev event=step_start step=score_$KIND
    # env -u LD_LIBRARY_PATH: the CUDA wheel paths set for Vipsania above break the scorer env (scipy import failed, 2026-11-06)
    env -u LD_LIBRARY_PATH CUDA_VISIBLE_DEVICES="" nice -n 10 taskset -c $SCPUS $PYG -I $HERE/vip_score.py $KIND --out $M/results/$KIND "${ftarg[@]}" --pre "${pre[@]}" >> $LOGS/m38-vipsania-score_$KIND.out 2>> $LOGS/m38-vipsania-score_$KIND.err
    local rc=$?
    ev event=step_end step=score_$KIND rc=$rc
    return $rc
  }
  score dev $DEV_GENOMES || { ev event=devb_aborted stage=score_dev; exit 23; }
  score setb $B_GENOMES  || { ev event=devb_aborted stage=score_setb; exit 24; }
fi
python3.12 -I - <<PY > $M/DEVB_DONE
import json, os, time
M = "$M"
out = dict(time=time.strftime("%Y-%m-%dT%H:%M:%S%z"), model="$MODEL_ID", runs=M + "/runs", results=M + "/results", events=M + "/events.jsonl")
for kind in ("dev", "setb"):
    p = f"{M}/results/{kind}/vip_{kind}_result.json"
    if os.path.exists(p):
        r = json.load(open(p))
        out[kind] = {n: (v["pooled"]["F1"] if kind == "dev" else v["table"]["pooled"]["F1"]) for n, v in r["variants"].items()}
print(json.dumps(out))
PY
ev event=devb_done
exit 0
