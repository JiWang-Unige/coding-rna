#!/bin/bash
# Stage driver of one M40 arm = m34_stage.sh with the M40 trainer (resumable; ARM_DONE marks the end).  Everything after training is the UNCHANGED M34 path:
#   stops 768 / 1536 / 2304 / 3072 -> m34_infer_retry.sh + m34_val.py [select: Ath NC_003071.7 + Osa NC_089044.1; own_species: the 8 new-species pseudo-contigs]
#   -> m34_select.py select (argmax Ath/Osa pooled F1 over {1536, 2304, 3072}, ties -> larger step) -> DEV exactly once (m34_select.py dev_start sentinel,
#   M32 c1_eval / c1_boot).  No transfer set (all 10 species are in every M40 arm).
# usage: m40_stage.sh <arm>     output: ~/coding-rna/outputs/M34-NTV2-SCALEUP/pc2_<arm>/
ARM=$1
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
ROOT=$HOME/coding-rna-reorg/experiments
O=$OUTROOT/$ARMPFX$ARM
mkdir -p $O $HOME/coding-rna/logs
exec 9>$O/.stage.lock
flock -n 9 || { echo "stage $ARM already running"; exit 3; }
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONUNBUFFERED=1
ev() { echo "{\"event\":\"$1\",\"arm\":\"$ARM\",\"time\":\"$(date -Is)\"$2}" >> $O/stage_log.jsonl; }
SEL=NC_003071.7,NC_089044.1
OWN=$(cd $M40S && $PY -c "import m40_common as W; print(','.join(W.T.arm_val_regions('$ARM')[1]))" 2>/dev/null)
[ -n "$OWN" ] || { ev own_regions_failed; exit 9; }
ev stage_start
cd $M34S
if [ ! -f $O/ARM_DONE ]; then
  for STOP in 768 1536 2304 3072; do
    if [ ! -f $O/DONE_$STOP.json ]; then
      bash $M40S/m40_supervise.sh $ARM $O $STOP || { ev train_failed ",\"stop\":$STOP"; exit 9; }
    fi
    if [ ! -f $O/valdone_$STOP ]; then
      V=$O/val_$STOP
      bash m34_infer_retry.sh $V $O/model_$(printf %06d $STOP).pt $SEL,$OWN 4 || { ev infer_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m34_val.py $V $STOP $SEL --tag select --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev val_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m34_val.py $V $STOP $OWN --tag own_species --arm $ARM --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || ev own_val_failed ",\"stop\":$STOP"
      rm -rf $V/pot; touch $O/valdone_$STOP
      ev validated ",\"stop\":$STOP"
    fi
  done
  BEST=$($PY m34_select.py select $ARM) || { ev select_failed; exit 9; }
  ev selected ",\"best_step\":$BEST"
  MODEL=$O/model_$(printf %06d $BEST).pt
  if [ ! -f $O/DEV_DONE ]; then
    $PY m34_select.py dev_start $ARM; rc=$?
    if [ $rc -eq 11 ]; then ev dev_refused; exit 11; fi
    if [ $rc -eq 0 ]; then
      ev dev_start ",\"best_step\":$BEST"
      D=$O/devM34_$ARM
      bash m34_infer_retry.sh $D $MODEL NC_003074.8,NC_089041.1 4 || { ev dev_infer_failed; exit 9; }
      ( cd $ROOT/$ENC_EXP/scripts && taskset -c 18-23 $PY c1_eval.py $D M40pc2_${ARM}_step$BEST DEV --compare-dev >> $D/eval.log 2>&1 && \
        taskset -c 18-23 $PY c1_boot.py $D $D/bootstrap_dev.json >> $D/eval.log 2>&1 ) || { ev dev_eval_failed; exit 9; }
      rm -rf $D/pot
      touch $O/DEV_DONE
      ev dev_done
    fi
  fi
  touch $O/ARM_DONE
fi
ev arm_done
