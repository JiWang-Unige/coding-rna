#!/bin/bash
# Stage driver of one M39 arm (copy of m34_stage.sh; protocol 29 §4.4-4.7; single instance per arm; fully resumable -- re-running continues where it stopped):
#   train to 768 / 1536 / 2304 / 3072 (3,072-step cosine schedule), after each stop run inference + internal validation
#     [tag select: old-species S chromosomes (Drosophila 2R + Tribolium chr 6) pooled; tag own_species: the 8 new-species validation crops]
#   -> select the final checkpoint among {1536, 2304, 3072} by the S pooled F1 only (m39_select.py)
#   -> old-species held-out chromosomes D (Drosophila X + Tribolium chr 7) EXACTLY ONCE for the final checkpoint (sentinels d_started.json / D_DONE; tag heldout_d)
#   No plant data (DEV / sealed / panel) is read.  P4-I is computed afterwards from metrics.jsonl (own_species lines of the two selected checkpoints).
# usage: [M39_BACKEND=pc2] m39_stage.sh <arm>   (arms: I10_k50U, I10_k50F, I10_T)   output: ~/coding-rna/outputs/M39-INSECT-C2/<arm>/ (pc2: pc2_<arm>/) (ARM_DONE marks the end)
ARM=$1
ROOT=$HOME/coding-rna-reorg/experiments
cd $ROOT/M39-INSECT-C2/scripts
source ./m39_env.sh          # M39_BACKEND (ntv2 | pc2) -> PY, ENC_EXP, ARMPFX, OUTROOT
O=$OUTROOT/$ARMPFX$ARM
mkdir -p $O $HOME/coding-rna/logs
exec 9>$O/.stage.lock
flock -n 9 || { echo "stage $ARM already running"; exit 3; }
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONUNBUFFERED=1
ev() { echo "{\"event\":\"$1\",\"arm\":\"$ARM\",\"time\":\"$(date -Is)\"$2}" >> $O/stage_log.jsonl; }
SEL=$($PY -c "import m39_train_common as T; print(','.join(T.select_regions()))" 2>/dev/null)
OWN=$($PY -c "import m39_train_common as T; print(','.join(T.arm_val_regions('$ARM')[1]))" 2>/dev/null)
DREG=$($PY -c "import m39_train_common as T; print(','.join(T.d_regions()))" 2>/dev/null)
[ -n "$SEL" ] && [ -n "$OWN" ] && [ -n "$DREG" ] || { ev regions_missing; exit 9; }
ev stage_start ",\"select\":\"$SEL\",\"own\":\"$OWN\",\"d\":\"$DREG\""
if [ ! -f $O/ARM_DONE ]; then
  for STOP in 768 1536 2304 3072; do
    if [ ! -f $O/DONE_$STOP.json ]; then
      bash m39_supervise.sh $ARM $O $STOP || { ev train_failed ",\"stop\":$STOP"; exit 9; }
    fi
    if [ ! -f $O/valdone_$STOP ]; then
      V=$O/val_$STOP
      bash m39_infer_retry.sh $V $O/model_$(printf %06d $STOP).pt $SEL,$OWN 4 || { ev infer_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m39_val.py $V $STOP $SEL --tag select --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev val_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m39_val.py $V $STOP $OWN --tag own_species --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev own_val_failed ",\"stop\":$STOP"; exit 9; }
      rm -rf $V/pot; touch $O/valdone_$STOP
      ev validated ",\"stop\":$STOP"
    fi
  done
  BEST=$($PY m39_select.py select $ARM) || { ev select_failed; exit 9; }
  ev selected ",\"best_step\":$BEST"
  MODEL=$O/model_$(printf %06d $BEST).pt
  if [ ! -f $O/D_DONE ]; then
    $PY m39_select.py d_start $ARM; rc=$?
    if [ $rc -eq 11 ]; then ev d_refused; exit 11; fi
    if [ $rc -eq 0 ]; then
      ev d_start ",\"best_step\":$BEST"
      D=$O/heldout_d
      bash m39_infer_retry.sh $D $MODEL $DREG 4 || { ev d_infer_failed; exit 9; }
      taskset -c 18-23 $PY m39_val.py $D $BEST $DREG --tag heldout_d --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev d_val_failed; exit 9; }
      rm -rf $D/pot
      touch $O/D_DONE
      ev d_done
    fi
  fi
  touch $O/ARM_DONE
fi
ev arm_done
