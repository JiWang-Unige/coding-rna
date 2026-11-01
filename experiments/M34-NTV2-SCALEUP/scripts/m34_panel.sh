#!/bin/bash
# P5 independent-panel inference driver (runs on ont, detached; single instance; every step resumable).  Outline paper_outline_20261027 P5; scope = core 10 species
# (report 25 §9); user decision gpu-after-k50f (2026-11-01): runs after M34 k50F, before the +2 seeds of T2 and k25U.
#   model: the final model T2 (pc2) at its selected checkpoint (2304, registered); recorded on first use, a different file is refused later.
#   per species (smallest first): m34_panel_infer.py over all its main sequences (zero-shot, same inference / decoding path as M34 and P2b, no posterior filter in the
#   GFF), potentials kept (--keep-head) for the CPU posterior queue m34_panel_postq.py (tau 1, as M33 / M36), which deletes each species' pot/ once its posteriors exist.
#   Nothing here reads a reference annotation and nothing is scored (scoring comes with the P6 baselines).
# Starts nothing on the GPU while a coding-rna GPU process or another M34 stage / inference of this user is alive (3 consecutive OK checks, 60 s apart).
# Events: $P/panel.jsonl ; sentinel: $P/PANEL_INFER_DONE ; outputs: $P/T2/<species>/<contig>/prediction/full_grammar.gff3, $P/post_T2/<species>/tau1/*.json
# usage: cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && bash m34_panel.sh       (normally started by m34_after_k50f.sh)
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
P=$O/panel
LOGS=$HOME/coding-rna/logs
mkdir -p $P $LOGS
exec 8>$P/.panel.lock
flock -n 8 || { echo "panel already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
export M34_BACKEND=pc2
source ./m34_env.sh                                   # -> PY = ~/envs/plantcad2/bin/python (pc2), OUTROOT
source $HOME/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts/m36_gate.sh      # m36_gpu_blockers
PYG=$HOME/envs/generanno/bin/python                   # CPU posterior queue (M33 / M36 environment)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
PCPUS=${PANEL_POST_CPUS:-24-31}
MIN_FREE_GB=${PANEL_MIN_FREE_GB:-250}
ME=$(id -un)
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $P/panel.jsonl; }
SPECIES=$($PY -c "import m34_panel_common as PC; print(' '.join(sorted(PC.CORE10, key=PC.total_bp)))") || { ev refused ",\"reason\":\"panel registry unreadable\""; exit 1; }

# ---- gate: no other coding-rna GPU job ---------------------------------------------------------------------------------------------------------------------
blockers() {
  local r="" pp g
  pp=$(pgrep -u $ME -f 'm34_trai[n]|m34_stag[e]|m34_supervis[e]|m34_p2[b]|m3[6]_queue|m34_infer\.p[y]' | tr '\n' ' ')
  [ -n "$pp" ] && r="$r processes: $pp;"
  g=$(m36_gpu_blockers); [ -n "$g" ] && r="$r gpu pids: $g;"
  echo "$r" | tr -d '"\\'
}
ev queued ",\"species\":\"$SPECIES\""
ok=0; n=0
while [ $ok -lt 3 ]; do
  R=$(blockers)
  if [ -z "$R" ]; then ok=$((ok+1)); sleep 60
  else ok=0; [ $((n % 60)) -eq 0 ] && ev waiting ",\"reasons\":\"$R\""; n=$((n+1)); sleep 300; fi
done

# ---- model: T2 at its registered selected checkpoint ------------------------------------------------------------------------------------------------------
B2=$(head -c 20 $O/pc2_T2/best_step.txt | tr -dc 0-9)
[ "$B2" = 2304 ] || { ev refused ",\"reason\":\"T2 selected step is $B2, registered 2304\""; exit 1; }
M2=$O/pc2_T2/model_$(printf %06d $B2).pt
[ -f $M2 ] || { ev refused ",\"reason\":\"model file missing\""; exit 1; }
h=$(md5sum $M2 | cut -d' ' -f1)
if [ -f $P/model_T2.json ]; then
  grep -q "\"md5\":\"$h\"" $P/model_T2.json || { ev refused ",\"reason\":\"T2 model changed since first use\""; exit 11; }
else
  echo "{\"arm\":\"T2\",\"model\":\"$M2\",\"md5\":\"$h\",\"bytes\":$(stat -c %s $M2),\"time\":\"$(date -Is)\"}" > $P/model_T2.json
fi
ev gate_open ",\"model\":\"$M2\",\"step\":$B2"

# ---- CPU posterior queue (started once; restarted if it exited on its idle timeout) -------------------------------------------------------------------------
postq() {
  pgrep -u $ME -f 'm34_panel_post[q]' > /dev/null && return 0
  setsid nohup nice -n 10 $PYG m34_panel_postq.py --run-root $P/T2 --post-root $P/post_T2 --workers 4 --cpus $PCPUS --mem-gb 30 --idle-timeout-h 48 \
    --log $P/postq.jsonl >> $LOGS/m34-panel-postq.out 2>&1 8>&- < /dev/null &
  ev postq_started ",\"pid\":$!"
}

# ---- inference per species (resumable per (contig, strand); retries like m34_infer_retry.sh) ----------------------------------------------------------------
for SP in $SPECIES; do
  OUT=$P/T2/$SP
  mkdir -p $OUT
  if ls $OUT/DONE_*.json >/dev/null 2>&1; then ev species_already_done ",\"species\":\"$SP\""; continue; fi
  while :; do       # disk guard: the posterior queue deletes pot/ species by species
    fr=$(df -BG --output=avail $P | tail -1 | tr -dc 0-9)
    [ "${fr:-0}" -ge $MIN_FREE_GB ] && break
    ev disk_wait ",\"free_gb\":$fr,\"need_gb\":$MIN_FREE_GB"; postq; sleep 1800
  done
  postq
  REG=$($PY -c "import m34_panel_common as PC; print(','.join(PC.contigs('$SP')))")
  ev infer_start ",\"species\":\"$SP\",\"contigs\":$(echo $REG | tr ',' '\n' | wc -l),\"bp\":$($PY -c "import m34_panel_common as PC; print(PC.total_bp('$SP'))")"
  done_sp=0
  for a in $(seq 1 60); do
    free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    if [ "${free:-0}" -lt 5000 ]; then sleep 120; continue; fi
    taskset -c $CPUS $PY m34_panel_infer.py --model $M2 --regions $REG --out $OUT --decode-workers 4 --cpus $CPUS --keep-head \
      >> $LOGS/m34-panel-$SP.out 2>> $LOGS/m34-panel-$SP.err < /dev/null
    rc=$?
    if ls $OUT/DONE_*.json >/dev/null 2>&1; then done_sp=1; break; fi
    ev infer_interruption ",\"species\":\"$SP\",\"attempt\":$a,\"rc\":$rc,\"stderr_tail\":\"$(tail -n 2 $LOGS/m34-panel-$SP.err | tr '\n"\\' '   ' | cut -c1-300)\""
    sleep 30
  done
  [ $done_sp = 1 ] || { ev infer_failed ",\"species\":\"$SP\""; exit 9; }
  ev infer_done ",\"species\":\"$SP\",\"summary\":$(cat $(ls $OUT/DONE_*.json | head -1))"
done
postq
ev panel_infer_done
touch $P/PANEL_INFER_DONE
