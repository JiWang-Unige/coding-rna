#!/bin/bash
# Stage driver of one M34 arm (single instance per arm; fully resumable -- re-running continues where it stopped):
#   train to 768 / 1536 / 2304 / 3072 (3,072-step cosine schedule), after each stop run inference + internal validation
#     [tag select: Ath NC_003071.7 + Osa NC_089044.1 pooled; tag own_species: the arm's new-species validation pseudo-contigs]
#   -> select the final checkpoint among {1536, 2304, 3072} by the Ath/Osa pooled F1 only (m34_select.py)
#   -> cross-species transfer set: validation pseudo-contigs of the species NOT in the arm (final checkpoint, once)
#   -> DEV (NC_003074.8 + NC_089041.1) EXACTLY ONCE for the final checkpoint (sentinels dev_started.json / DEV_DONE), eval + 1-Mb block bootstrap (c1_eval / c1_boot of the
#      backend's experiment dir, M31 or M32, unchanged).
# usage: [M34_BACKEND=pc2] m34_stage.sh <arm>   output: ~/coding-rna/outputs/M34-NTV2-SCALEUP/<arm>/ (pc2: pc2_<arm>/)  (ARM_DONE marks the end)
ARM=$1
ROOT=$HOME/coding-rna-reorg/experiments
cd $ROOT/M34-NTV2-SCALEUP/scripts
source ./m34_env.sh          # M34_BACKEND (ntv2 | pc2) -> PY, ENC_EXP, ARMPFX
O=$OUTROOT/$ARMPFX$ARM
mkdir -p $O $HOME/coding-rna/logs
exec 9>$O/.stage.lock
flock -n 9 || { echo "stage $ARM already running"; exit 3; }
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONUNBUFFERED=1
ev() { echo "{\"event\":\"$1\",\"arm\":\"$ARM\",\"time\":\"$(date -Is)\"$2}" >> $O/stage_log.jsonl; }
SEL=NC_003071.7,NC_089044.1
OWN=$($PY -c "import m34_train_common as T; print(','.join(T.arm_val_regions('$ARM')[1]))" 2>/dev/null)
OTHER=$($PY -c "import m34_train_common as T; print(','.join(T.other_val_regions('$ARM')))" 2>/dev/null)
ev stage_start
if [ ! -f $O/ARM_DONE ]; then
  for STOP in 768 1536 2304 3072; do
    if [ ! -f $O/DONE_$STOP.json ]; then
      bash m34_supervise.sh $ARM $O $STOP || { ev train_failed ",\"stop\":$STOP"; exit 9; }
    fi
    if [ ! -f $O/valdone_$STOP ]; then
      V=$O/val_$STOP
      REG=$SEL; [ -n "$OWN" ] && REG=$SEL,$OWN
      bash m34_infer_retry.sh $V $O/model_$(printf %06d $STOP).pt $REG 4 || { ev infer_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m34_val.py $V $STOP $SEL --tag select --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev val_failed ",\"stop\":$STOP"; exit 9; }
      if [ -n "$OWN" ]; then taskset -c 18-23 $PY m34_val.py $V $STOP $OWN --tag own_species --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || ev own_val_failed ",\"stop\":$STOP"; fi
      rm -rf $V/pot; touch $O/valdone_$STOP
      ev validated ",\"stop\":$STOP"
    fi
  done
  BEST=$($PY m34_select.py select $ARM) || { ev select_failed; exit 9; }
  ev selected ",\"best_step\":$BEST"
  MODEL=$O/model_$(printf %06d $BEST).pt
  if [ -n "$OTHER" ] && [ ! -f $O/valdone_transfer ]; then
    bash m34_infer_retry.sh $O/transfer $MODEL $OTHER 4 || { ev transfer_infer_failed; exit 9; }
    taskset -c 18-23 $PY m34_val.py $O/transfer $BEST $OTHER --tag transfer --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || ev transfer_val_failed
    rm -rf $O/transfer/pot; touch $O/valdone_transfer
    ev transfer_done
  fi
  if [ ! -f $O/DEV_DONE ]; then
    $PY m34_select.py dev_start $ARM; rc=$?
    if [ $rc -eq 11 ]; then ev dev_refused; exit 11; fi
    if [ $rc -eq 0 ]; then
      ev dev_start ",\"best_step\":$BEST"
      D=$O/devM34_$ARM
      bash m34_infer_retry.sh $D $MODEL NC_003074.8,NC_089041.1 4 || { ev dev_infer_failed; exit 9; }
      ( cd $ROOT/$ENC_EXP/scripts && taskset -c 18-23 $PY c1_eval.py $D M34${M34_BACKEND}_${ARM}_step$BEST DEV --compare-dev >> $D/eval.log 2>&1 && \
        taskset -c 18-23 $PY c1_boot.py $D $D/bootstrap_dev.json >> $D/eval.log 2>&1 ) || { ev dev_eval_failed; exit 9; }
      rm -rf $D/pot
      touch $O/DEV_DONE
      ev dev_done
    fi
  fi
  touch $O/ARM_DONE
fi
ev arm_done
