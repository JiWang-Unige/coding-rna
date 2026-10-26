#!/bin/bash
# M36 GPU queue (runs on ont, detached; single instance; every stage resumable -- re-running continues where it stopped).  Protocol: docs/review_20261016/22_M36_bootstrap_adaptation_protocol.md
# (base model = M34 T10 step 3072, section 9); implementation: docs/review_20261016/26_M36_cpu_preparation.md section 7.
#   gate    : waits for  $M34O/p2b/P2B_DONE  (P2b B-set inference + scoring of T10 / T2 finished), no other coding-rna GPU job, 3 consecutive OK checks 60 s apart (P2b gate logic,
#             m36_gate.sh).  The TE project's GPU processes (same account) do not count.  Before every later GPU stage it waits again until no coding-rna GPU process is left.
#   A zero-shot : the base model on the label-set A contigs (train chromosomes + the two __val20M pseudo-contigs), pot/ kept            -> $O/zeroshot_A
#   CPU (detached, taskset 24-31, nice 10; starts with the zero-shot run and works on each contig as soon as its pot/chains exist):
#               m33 posterior (tau 1) + no-gene marginal per (contig, strand) -> S1 seeds (log P >= -1.400, CDS >= 96), P>=0.9 set, pseudo-background (P(no gene) >= 0.99)
#                                                                              -> $O/post_A, $O/nogene_A, $O/seeds/{seeds_S1_*, background_*, READY_<species>}
#   trainings : 7 x 512 steps, serially (S2 first: it needs no S1 seeds; an arm whose seeds are not ready yet is skipped for the next ready one)   -> $O/train/<arm>_<species>
#   B inference: one-shot per trained arm (model_000512.pt), sentinel $O/evalB/<arm>_<species>.B_INFER_STARTED.json                           -> $O/evalB/<arm>_<species>
#   scoring   : m36_score.py main (one-shot, sentinel results/SCORE_MAIN_STARTED.json) -> $O/results/m36_result.json, then M36_DONE.  Then (descriptive, CPU) the posterior-filtered
#               scoring m36_score.py post once the B posteriors ($O/postB) are done -> M36_POST_DONE.
# Z arm: NOT re-run -- the P2b T10 B-set output $M34O/p2b/T10 (same model, same path, pot/ kept).
# Events (one JSON line each): $O/queue.jsonl.  Final sentinel: $O/M36_DONE.  Exit codes: 3 already running, 9 a stage failed, 11 refused (inputs changed).
# usage: cd ~/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts && setsid nohup bash m36_queue.sh > ~/coding-rna/logs/m36-queue.out 2>&1 < /dev/null &
# test knobs: M36_OUT / M36_M34_OUT / M36_LOGS (fake roots), M36_DRY=1 (walk the queue logic, run nothing), M36_GATE_SKIP_PROCS=1, M36_SLEEP_OK / M36_SLEEP_WAIT, M36_ORDER, M36_POST_B=0
EXP=${M36_EXP:-$HOME/coding-rna-reorg/experiments}
S=$EXP/M36-BOOTSTRAP-ADAPT/scripts
M34S=$EXP/M34-NTV2-SCALEUP/scripts
O=${M36_OUT:-$HOME/coding-rna/outputs/M36-BOOTSTRAP-ADAPT}
M34O=${M36_M34_OUT:-$HOME/coding-rna/outputs/M34-NTV2-SCALEUP}
LOGS=${M36_LOGS:-$HOME/coding-rna/logs}
export M36_OUT=$O M36_M34_OUT=$M34O M36_LOGS=$LOGS M36_SCRIPTS=$S M36_M34_SCRIPTS=$M34S
mkdir -p $O $LOGS
exec 8>$O/.queue.lock
flock -n 8 || { echo "m36 queue already running"; exit 3; }
ME=$(id -un)
PYG=${M36_PY_CPU:-$HOME/envs/generanno/bin/python}
BASE=$M34O/pc2_T10/model_003072.pt
Z=$M34O/p2b/T10
P2B_DONE=$M34O/p2b/P2B_DONE
SLEEP_OK=${M36_SLEEP_OK:-60}; SLEEP_WAIT=${M36_SLEEP_WAIT:-300}
STEPS=512
# regions (M30 protocol 13 section 2): label set A = training chromosomes + validation pseudo-contig; evaluation set B per species
A_REGIONS=NC_141366.1,NC_141368.1,NC_141370.1,NC_141372.1,NC_141374.1__val20M,NW_028238198.1,NW_028238420.1,NW_028238472.1__val20M
declare -A B_REGIONS=( [dryas_octopetala]=NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1 [eschscholzia_californica]=NW_028238309.1,NW_028238461.1,NW_028238483.1 )
source $S/m36_gate.sh
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $O/queue.jsonl; }
die() { ev "$1" "$2"; exit ${3:-9}; }
dry() { [ "${M36_DRY:-0}" = 1 ]; }

# ---- gate -----------------------------------------------------------------------------------------------------------------------------------------
wait_gate() {
  local ok=0 n=0 R
  while [ $ok -lt 3 ]; do
    R=$(m36_blockers)
    if [ -z "$R" ]; then ok=$((ok+1)); sleep $SLEEP_OK
    else
      ok=0
      [ $((n % 60)) -eq 0 ] && ev waiting ",\"reasons\":\"$R\""       # first time and then every 5 h
      n=$((n+1)); sleep $SLEEP_WAIT
    fi
  done
}
gpu_clear_wait() {            # before every GPU stage: no coding-rna GPU process of this user
  [ "${M36_GATE_SKIP_PROCS:-0}" = 1 ] && return 0
  local g n=0
  while :; do
    g=$(m36_gpu_blockers)
    [ -z "$g" ] && return 0
    [ $((n % 30)) -eq 0 ] && ev gpu_wait ",\"pids\":\"$g\""
    n=$((n+1)); sleep 120
  done
}

# ---- base model / Z arm ---------------------------------------------------------------------------------------------------------------------------
verify_base() {
  [ -f $BASE ] || die base_missing ",\"path\":\"$BASE\""
  local h; h=$(md5sum $BASE | cut -d' ' -f1)
  grep -Eq "\"md5\": ?\"$h\"" $M34O/p2b/model_T10.json || die base_md5_mismatch ",\"md5\":\"$h\",\"reason\":\"model differs from the one P2b used for T10 (= the Z arm)\"" 11
  [ "$(head -c 20 $M34O/pc2_T10/best_step.txt | tr -dc 0-9)" = 3072 ] || die base_step_mismatch ",\"reason\":\"pc2_T10 selected step is not 3072\"" 11
  ls $Z/DONE_*.json >/dev/null 2>&1 || die z_arm_missing ",\"reason\":\"no DONE_*.json in $Z\""
  [ -d $Z/pot ] || ev z_pot_missing ",\"note\":\"Z posterior-filtered descriptive scoring will not be possible\""
  if [ -f $O/base_model.json ]; then
    grep -Eq "\"md5\": ?\"$h\"" $O/base_model.json || die base_changed ",\"md5\":\"$h\"" 11
  else
    echo "{\"model\":\"$BASE\",\"md5\":\"$h\",\"bytes\":$(stat -c %s $BASE),\"z_dir\":\"$Z\",\"time\":\"$(date -Is)\"}" > $O/base_model.json
  fi
  ev base_verified ",\"md5\":\"$h\""
}

# ---- CPU queues (detached) ------------------------------------------------------------------------------------------------------------------------
cpuq_alive() { pgrep -u $ME -f "m36_cpu_queue.p[y] $1" >/dev/null; }
SEEDQ_STARTS=0
start_cpuq_seeds() {
  if dry; then ev would_run ",\"stage\":\"cpu_seeds_queue\""; return 0; fi
  cpuq_alive seeds && return 0
  local sp all=1
  for sp in dryas_octopetala eschscholzia_californica; do [ -f $O/seeds/READY_$sp ] || all=0; done
  [ $all = 1 ] && return 0
  SEEDQ_STARTS=$((SEEDQ_STARTS+1))
  [ $SEEDQ_STARTS -le 6 ] || die cpu_seeds_failed ",\"reason\":\"CPU seed queue restarted 6 times without producing the S1 seeds; see $O/cpuq_seeds.jsonl\""
  setsid nohup nice -n 10 $PYG $S/m36_cpu_queue.py seeds --run $O/zeroshot_A --post $O/post_A --nogene $O/nogene_A --seeds-out $O/seeds --workers 8 --cpus 24-31 --mem-gb 45 \
    --idle-timeout-h 48 --log $O/cpuq_seeds.jsonl >> $LOGS/m36-cpuq-seeds.out 2>&1 8>&- < /dev/null &
  ev cpu_seeds_started ",\"pid\":$!"
}
start_cpuq_bpost() {
  [ "${M36_POST_B:-1}" = 1 ] || return 0
  if dry; then ev would_run ",\"stage\":\"cpu_bpost_queue\""; return 0; fi
  cpuq_alive bpost && return 0
  [ -f $O/M36_POST_DONE ] && return 0
  local runs=() arm sp
  for sp in dryas_octopetala eschscholzia_californica; do runs+=("Z_$sp=$Z:$sp:keep"); done
  while read -r arm sp; do runs+=("${arm}_${sp}=$O/evalB/${arm}_${sp}:$sp"); done < <(order_lines)
  setsid nohup nice -n 10 $PYG $S/m36_cpu_queue.py bpost --runs "${runs[@]}" --post-root $O/postB --workers 4 --cpus 24-31 --mem-gb 30 --idle-timeout-h 72 \
    --log $O/cpuq_bpost.jsonl >> $LOGS/m36-cpuq-bpost.out 2>&1 8>&- < /dev/null &
  ev cpu_bpost_started ",\"pid\":$!"
}
order_lines() { ( cd $S && $PYG -c "import m36_common as W; print('\n'.join(a + ' ' + s for a, s in W.order()))" ); }

# ---- stages ---------------------------------------------------------------------------------------------------------------------------------------
zeroshot_A() {
  [ -f $O/A_ZEROSHOT_DONE ] && return 0
  if dry; then ev would_run ",\"stage\":\"zeroshot_A\",\"regions\":\"$A_REGIONS\""; return 0; fi
  gpu_clear_wait
  ev zeroshot_start ",\"regions\":\"$A_REGIONS\""
  M36_KEEP_POT=1 bash $S/m36_infer.sh $O/zeroshot_A $BASE $A_REGIONS 4 || die zeroshot_failed
  touch $O/A_ZEROSHOT_DONE
  ev zeroshot_done ",\"dir\":\"$O/zeroshot_A\""
}
arm_ready() {                 # S2 needs only the protein seeds; S1 / S3 / F1 need the S1 seeds + background of the species
  case $1 in S2) return 0;; *) [ -f $O/seeds/READY_$2 ];; esac
}
train_arm() {
  local arm=$1 sp=$2 key=${1}_${2}; local D=$O/train/$key
  [ -f $D/DONE_$STEPS.json ] && return 0
  if dry; then ev would_run ",\"stage\":\"train\",\"arm\":\"$key\",\"steps\":$STEPS"; return 0; fi
  gpu_clear_wait
  ev train_start ",\"arm\":\"$key\""
  taskset -c 24-31 nice -n 10 $PYG $S/m36_schedule.py make $arm $sp >> $LOGS/m36-schedule-$key.out 2>&1 || die schedule_failed ",\"arm\":\"$key\""
  bash $S/m36_supervise.sh $arm $sp $D $STEPS; local rc=$?
  [ $rc -eq 0 ] || die train_failed ",\"arm\":\"$key\",\"rc\":$rc" $rc
  ev train_done ",\"arm\":\"$key\""
}
infer_b() {
  local arm=$1 sp=$2 key=${1}_${2}; local model=$O/train/$key/model_$(printf %06d $STEPS).pt run=$O/evalB/$key
  ls $run/DONE_*.json >/dev/null 2>&1 && return 0
  if dry; then ev would_run ",\"stage\":\"infer_B\",\"arm\":\"$key\",\"regions\":\"${B_REGIONS[$sp]}\""; return 0; fi
  [ -f $model ] || die model_missing ",\"arm\":\"$key\""
  local h; h=$(md5sum $model | cut -d' ' -f1)
  $PYG $S/m36_sentinel.py start $O/evalB/$key.B_INFER_STARTED.json "{\"model\":\"$model\",\"md5\":\"$h\",\"regions\":\"${B_REGIONS[$sp]}\"}"; local rc=$?
  [ $rc -eq 0 ] || die b_infer_refused ",\"arm\":\"$key\",\"rc\":$rc" 11
  gpu_clear_wait
  ev infer_B_start ",\"arm\":\"$key\",\"model_md5\":\"$h\""
  M36_KEEP_POT=${M36_POST_B:-1} bash $S/m36_infer.sh $run $model ${B_REGIONS[$sp]} 4 || die infer_B_failed ",\"arm\":\"$key\""
  $PYG $S/m36_sentinel.py done $O/evalB/$key.B_INFER_STARTED.json
  ev infer_B_done ",\"arm\":\"$key\""
}
score_main() {
  if [ -f $O/M36_DONE ]; then return 0; fi
  if dry; then ev would_run ",\"stage\":\"score_main\""; return 0; fi
  ev score_start
  taskset -c 18-23 nice -n 10 $PYG $S/m36_score.py main --root $O >> $LOGS/m36-score.out 2>&1; local rc=$?
  if [ $rc -ne 0 ] && { [ $rc -ne 10 ] || [ ! -f $O/results/m36_result.json ]; }; then die score_failed ",\"rc\":$rc" $rc; fi
  local SUM
  SUM=$($PYG -c "
import json, sys
r = json.load(open(sys.argv[1])); v = r['verdict']
print(json.dumps(dict(label=v['label'], go_strong=v['go_strong'], go=v['go'], no_go=v['no_go'], c2_retracted=v['c2_retracted'],
    F1_Z={s: round(x['F1_Z'], 2) for s, x in r['species'].items()}, F1_CarbonA_filt={s: round(x['F1_CarbonA_filt'], 2) for s, x in r['species'].items()},
    arms={s: {a: dict(F1=round(y['F1'], 2), closure=round(y['closure'], 3)) for a, y in x['arms'].items()} for s, x in r['species'].items()}, result=sys.argv[1])))" $O/results/m36_result.json) || die summary_failed
  touch $O/M36_DONE
  ev m36_done ",\"summary\":$SUM"
}
score_post() {
  [ "${M36_POST_B:-1}" = 1 ] || return 0
  [ -f $O/M36_POST_DONE ] && return 0
  if dry; then ev would_run ",\"stage\":\"score_post\""; return 0; fi
  while cpuq_alive bpost; do sleep 600; done
  ev score_post_start
  taskset -c 18-23 nice -n 10 $PYG $S/m36_score.py post --root $O >> $LOGS/m36-score-post.out 2>&1; local rc=$?
  if [ $rc -eq 0 ] || [ $rc -eq 10 ]; then touch $O/M36_POST_DONE; ev m36_post_done; else ev score_post_failed ",\"rc\":$rc"; fi
}

# ---- run ------------------------------------------------------------------------------------------------------------------------------------------
ev queued ",\"dry\":${M36_DRY:-0}"
wait_gate
ev gate_open
verify_base
mapfile -t ORDER < <(order_lines | tr ' ' ':')
ev plan ",\"order\":\"${ORDER[*]}\",\"steps\":$STEPS"
start_cpuq_seeds
zeroshot_A
PENDING=("${ORDER[@]}")
nwait=0; VDONE=" "                                     # VDONE: arms 'trained' in a dry run
while [ ${#PENDING[@]} -gt 0 ]; do
  PICK=""; NEW=()
  for item in "${PENDING[@]}"; do
    arm=${item%%:*}; sp=${item#*:}
    if [ -f $O/train/${arm}_${sp}/DONE_$STEPS.json ]; then continue; fi
    case "$VDONE" in *" $item "*) continue;; esac
    NEW+=("$item")
    [ -z "$PICK" ] && { dry || arm_ready $arm $sp; } && PICK=$item
  done
  PENDING=("${NEW[@]}")
  [ ${#PENDING[@]} -eq 0 ] && break
  if [ -z "$PICK" ]; then
    [ $((nwait % 12)) -eq 0 ] && ev waiting_seeds ",\"pending\":\"${PENDING[*]}\""
    nwait=$((nwait+1))
    start_cpuq_seeds                                 # restarts the CPU queue if it died before the seeds were ready
    sleep $SLEEP_WAIT; continue
  fi
  train_arm ${PICK%%:*} ${PICK#*:}
  dry && VDONE="$VDONE$PICK "
done
ev trainings_done
start_cpuq_bpost
for item in "${ORDER[@]}"; do infer_b ${item%%:*} ${item#*:}; done
ev inference_B_done
score_main
score_post
ev queue_exit
