#!/bin/bash
# M34 P2b driver (runs on ont, detached; single instance; every step resumable).  Protocol: docs/review_20261016/19_M34_scaleup_plan.md section 3.1, bullet P2b:
#   zero-shot exact-chain F1 of T10 and T2 (pc2 backend, each with its selected checkpoint) on the M30 evaluation set B (7 contigs, 273 Mb), same inference/decoding path as M34
#   (m34_infer via m34_p2b_infer.py), no posterior filter, then T10 - T2 with the paired 1-Mb block bootstrap (m34_p2b_score.py).
# Starts NOTHING on the GPU until  $O/M34_STAGE1_DONE  exists AND no m34_train / m34_infer / m34_stage / m36 process of this user is left and no GPU compute process of this user
# exists (3 consecutive OK checks, 60 s apart).  Events (one JSON line each): $O/p2b/p2b.jsonl ; final sentinel: $O/p2b/P2B_DONE ; result: $O/p2b/p2b_result.json.
# The T10 B-set run directory ($O/p2b/T10) is the zero-shot arm Z of M36 (its potentials pot/ are kept: P2B_KEEP_POT_T10=1 default; T2's are deleted after scoring).
# usage: cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && setsid nohup bash m34_p2b.sh > ~/coding-rna/logs/m34-p2b.out 2>&1 < /dev/null &
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
P=$O/p2b
LOGS=$HOME/coding-rna/logs
mkdir -p $P $LOGS
exec 8>$P/.p2b.lock
flock -n 8 || { echo "p2b already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
export M34_BACKEND=pc2
source ./m34_env.sh                                   # -> PY = ~/envs/plantcad2/bin/python (pc2), OUTROOT
PYG=$HOME/envs/generanno/bin/python                   # CPU scoring: the M30 ruler's environment (as M33 m33_evalB)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
SCPUS=${P2B_SCORE_CPUS:-18-23}
ME=$(id -un)
REG=NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1,NW_028238309.1,NW_028238461.1,NW_028238483.1     # evaluation set B (13_M30_protocol section 2)
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $P/p2b.jsonl; }

# ---- gate: STAGE1 finished and no GPU job of this user -------------------------------------------------------------------------------------------------
blockers() {   # prints the reasons that block the start (empty = free).  The patterns contain a bracket so that pgrep never matches its own command line.
  local r=""
  [ -f $O/M34_STAGE1_DONE ] || r="no M34_STAGE1_DONE;"
  local pp; pp=$(pgrep -u $ME -f 'm34_trai[n]|m34_infe[r]|m34_stag[e]|m34_queu[e]|m34_supervis[e]|m34_smok[e]|m3[6]_|M3[6]-BOOTSTRAP' | tr '\n' ' ')
  [ -n "$pp" ] && r="$r processes: $pp;"
  local g=""   # GPU compute processes of this user that belong to coding-rna (the ab-initio-TE project runs under the same account and also uses the GPU: it does not count)
  for p in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    [ "$(ps -o user= -p $p 2>/dev/null | tr -d ' ')" = "$ME" ] || continue
    if tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -q 'coding-rna' || readlink /proc/$p/cwd 2>/dev/null | grep -q 'coding-rna'; then g="$g $p"; fi
  done
  [ -n "$g" ] && r="$r gpu pids:$g;"
  echo "$r" | tr -d '"\\'
}
ev queued
ok=0; n=0
while [ $ok -lt 3 ]; do
  R=$(blockers)
  if [ -z "$R" ]; then ok=$((ok+1)); sleep 60
  else
    ok=0
    [ $((n % 60)) -eq 0 ] && ev waiting ",\"reasons\":\"$R\""      # first time and then every 5 h
    n=$((n+1)); sleep 300
  fi
done

# ---- models (selected checkpoints; T10 is read at run time) --------------------------------------------------------------------------------------------
[ -f $O/pc2_T10/ARM_DONE ] && [ -f $O/pc2_T10/selection.json ] && [ -f $O/pc2_T10/best_step.txt ] || { ev refused ",\"reason\":\"pc2_T10 not finished (ARM_DONE / selection.json / best_step.txt)\""; exit 1; }
B10=$(head -c 20 $O/pc2_T10/best_step.txt | tr -dc 0-9); B2=$(head -c 20 $O/pc2_T2/best_step.txt | tr -dc 0-9)
[ "$B2" = 2304 ] || { ev refused ",\"reason\":\"T2 selected step is $B2, registered 2304\""; exit 1; }
M10=$O/pc2_T10/model_$(printf %06d $B10).pt; M2=$O/pc2_T2/model_$(printf %06d $B2).pt
[ -f $M10 ] && [ -f $M2 ] || { ev refused ",\"reason\":\"model file missing\""; exit 1; }
fp() {   # arm model : record the model on first use, refuse a different file later (the same model must serve every restart)
  local h; h=$(md5sum $2 | cut -d' ' -f1)
  if [ -f $P/model_$1.json ]; then
    grep -q "\"md5\":\"$h\"" $P/model_$1.json || { ev refused ",\"reason\":\"model of $1 changed since first use\""; exit 11; }
  else
    echo "{\"arm\":\"$1\",\"model\":\"$2\",\"md5\":\"$h\",\"bytes\":$(stat -c %s $2),\"time\":\"$(date -Is)\"}" > $P/model_$1.json
  fi
}
fp T10 $M10; fp T2 $M2
ev gate_open ",\"T10_step\":$B10,\"T10_model\":\"$M10\",\"T2_step\":$B2,\"T2_model\":\"$M2\""

# ---- inference (resumable per (contig, strand); retries like m34_infer_retry.sh) -------------------------------------------------------------------------
infer() {   # arm model outdir [extra m34_p2b_infer args]
  local ARM=$1 MODEL=$2 OUT=$3; shift 3
  mkdir -p $OUT
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then ev infer_already_done ",\"arm\":\"$ARM\""; return 0; fi
  ev infer_start ",\"arm\":\"$ARM\",\"model\":\"$MODEL\",\"out\":\"$OUT\",\"regions\":\"$REG\""
  local a free rc
  for a in $(seq 1 60); do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "${free:-0}" -lt 5000 ]; then sleep 120; continue; fi
    taskset -c $CPUS $PY m34_p2b_infer.py --model $MODEL --regions $REG --out $OUT --decode-workers 4 --cpus $CPUS "$@" >> $LOGS/m34-p2b-$ARM.out 2>> $LOGS/m34-p2b-$ARM.err < /dev/null
    rc=$?
    if ls $OUT/DONE_*.json >/dev/null 2>&1; then ev infer_done ",\"arm\":\"$ARM\",\"attempt\":$a"; return 0; fi
    ev infer_interruption ",\"arm\":\"$ARM\",\"attempt\":$a,\"rc\":$rc,\"stderr_tail\":\"$(tail -n 2 $LOGS/m34-p2b-$ARM.err | tr '\n"\\' '   ' | cut -c1-300)\""
    sleep 30
  done
  return 1
}
KEEP10=""; [ "${P2B_KEEP_POT_T10:-1}" = 1 ] && KEEP10="--keep-head"
infer T10 $M10 $P/T10 $KEEP10 || { ev infer_failed ",\"arm\":\"T10\""; exit 9; }
infer T2 $M2 $P/T2 || { ev infer_failed ",\"arm\":\"T2\""; exit 9; }

# ---- scoring + paired 1-Mb block bootstrap (CPU) -------------------------------------------------------------------------------------------------------
ev score_start
taskset -c $SCPUS nice -n 10 $PYG m34_p2b_score.py --new T10=$P/T10 --old T2=$P/T2 --out $P/p2b_result.json --label P2b --z-arm T10 >> $LOGS/m34-p2b-score.out 2>&1 \
  || { ev score_failed; exit 9; }
rm -rf $P/T2/pot
SUM=$($PYG -c "
import json, sys
r = json.load(open(sys.argv[1])); v = r['verdict']
print(json.dumps(dict(met=v['met'], checks=v['checks'], diff_point=v['point'], diff_CI95=v['CI95'],
                      F1={n: {sc: round(x['table'][sc]['F1'], 2) for sc in ('pooled', 'dryas_octopetala', 'eschscholzia_californica')} for n, x in r['runs'].items()},
                      T10_step=$B10, T2_step=$B2, z_arm_dir='$P/T10', result='$P/p2b_result.json')))" $P/p2b_result.json) || { ev summary_failed; exit 9; }
ev p2b_done ",\"summary\":$SUM"
touch $P/P2B_DONE
