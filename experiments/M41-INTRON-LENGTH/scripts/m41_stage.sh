#!/bin/bash
# M41 stage driver (protocol docs/review_20261016/33 section 3; = M34 m34_stage.sh for a legacy-only arm, with the M41 trainer / inference):
#   train T2_ILS to 768 / 1536 / 2304 / 3072 (T2 schedule file, seeds, 3,072-step cosine; resumable, crash-resilient as m34_supervise.sh)
#   -> after each stop: M41 inference + M34 internal validation (tag select: Ath NC_003071.7 + Osa NC_089044.1 pooled)
#   -> select among {1536, 2304, 3072} by the Ath/Osa pooled F1 only (m41_select.py = the M34 rule)
#   -> DEV (NC_003074.8 + NC_089041.1) EXACTLY ONCE (sentinels dev_started.json / DEV_DONE): c1_eval + c1_boot (M32, unchanged) + m41_compare (paired, by stratum)
#   -> evaluation set B zero-shot (m41_infer --setb) + m34_p2b_score (new T2_ILS vs old T2)  -> ARM_DONE
# usage: m41_stage.sh            (one instance; GPU queue placement is the maintenance session's job)
# env: M41_CPUS (default 2,4,6,8,10,12,14,16-24, as M34), M41_MIN_FREE_MIB (6000), M41_SHORT_STAGES (1, pre-registered)
LABEL=T2_ILS; ARM=T2
R=${M41_SHORT_STAGES:-1}
EXPS=$HOME/coding-rna-reorg/experiments
S41=$EXPS/M41-INTRON-LENGTH/scripts; S34=$EXPS/M34-NTV2-SCALEUP/scripts
O=$HOME/coding-rna/outputs/M41-INTRON-LENGTH/pc2_$LABEL
O34=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
LOGS=$HOME/coding-rna/logs
mkdir -p $O $LOGS
exec 9>$O/.stage.lock
flock -n 9 || { echo "stage $LABEL already running"; exit 3; }
export M34_BACKEND=pc2 M34_LORA=on HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1 \
       PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd $S34 && source ./m34_env.sh            # PY = plantcad2 env
GPY=$HOME/envs/generanno/bin/python       # ruler env (m34_p2b_score, m41_compare)
CPUS=${M41_CPUS:-2,4,6,8,10,12,14,16-24}
SEL=NC_003071.7,NC_089044.1
DEVR=NC_003074.8,NC_089041.1
REGB=NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1,NW_028238309.1,NW_028238461.1,NW_028238483.1
ev() { echo "{\"event\":\"$1\",\"arm\":\"$LABEL\",\"time\":\"$(date -Is)\"$2}" >> $O/stage_log.jsonl; }
free_mib() { nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1; }

train_to() {   # stop : resumable trainer with the m34_supervise.sh retry policy
  local STOP=$1 laststep=-1 noprog=0 attempt rc last err LOG=$LOGS/m41-train-$LABEL-$(date +%Y%m%d-%H%M%S)
  for attempt in $(seq 1 80); do
    if [ "$(free_mib)" -lt ${M41_MIN_FREE_MIB:-6000} ]; then echo "{\"event\":\"gpu_busy_wait\",\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl; sleep 120; continue; fi
    taskset -c $CPUS $PY $S41/m41_train.py run --arm $ARM --label $LABEL --out $O --steps 3072 --stop-at $STOP --batch 8 --workers 12 --cpus $CPUS \
      --short-stages $R --resume >> $LOG.out 2>> $LOG.err < /dev/null
    rc=$?
    [ -f $O/DONE_$STOP.json ] && return 0
    last=$($PY -c "import torch;print(torch.load('$O/ckpt_latest.pt',map_location='cpu',weights_only=False)['step'])" 2>/dev/null || echo -1)
    err=$(tail -n 3 $LOG.err | tr '\n"' '  ' | cut -c1-400)
    echo "{\"event\":\"interruption\",\"attempt\":$attempt,\"rc\":$rc,\"last_checkpoint_step\":\"$last\",\"stderr_tail\":\"$err\",\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl
    if [ "$last" = "$laststep" ]; then noprog=$((noprog+1)); else noprog=0; fi
    laststep=$last
    if [ $noprog -ge 3 ] && ! echo "$err" | grep -qiE 'out of memory|CUDA'; then echo "{\"event\":\"FATAL\",\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl; return 9; fi
    pkill -P $$ -f m41_dp 2>/dev/null; sleep 30
  done
  return 9
}

infer() {   # outdir model regions [--setb] : retrying, resumable per (region, strand)
  local D=$1 MODEL=$2 REG=$3 EXTRA=$4 a rc
  mkdir -p $D
  for a in $(seq 1 40); do
    if [ "$(free_mib)" -lt 5000 ]; then sleep 120; continue; fi
    taskset -c $CPUS $PY $S41/m41_infer.py $EXTRA --model $MODEL --regions $REG --out $D --decode-workers 4 --cpus $CPUS \
      >> $LOGS/m41-infer-$LABEL-$(basename $D).out 2>> $LOGS/m41-infer-$LABEL-$(basename $D).err < /dev/null
    rc=$?
    ls $D/DONE_*.json >/dev/null 2>&1 && return 0
    echo "{\"event\":\"infer_interruption\",\"dir\":\"$D\",\"attempt\":$a,\"rc\":$rc,\"time\":\"$(date -Is)\"}" >> $O/interruptions.jsonl
    sleep 30
  done
  return 9
}

ev stage_start ",\"short_stages\":$R"
if [ ! -f $O/ARM_DONE ]; then
  for STOP in 768 1536 2304 3072; do
    [ -f $O/DONE_$STOP.json ] || train_to $STOP || { ev train_failed ",\"stop\":$STOP"; exit 9; }
    if [ ! -f $O/valdone_$STOP ]; then
      V=$O/val_$STOP
      infer $V $O/model_$(printf %06d $STOP).pt $SEL || { ev infer_failed ",\"stop\":$STOP"; exit 9; }
      taskset -c 18-23 $PY m34_val.py $V $STOP $SEL --tag select --arm $LABEL --metrics $O/metrics.jsonl >> $O/val.log 2>&1 || { ev val_failed ",\"stop\":$STOP"; exit 9; }
      rm -rf $V/pot; touch $O/valdone_$STOP
      ev validated ",\"stop\":$STOP"
    fi
  done
  BEST=$($PY $S41/m41_select.py select $O $LABEL) || { ev select_failed; exit 9; }
  ev selected ",\"best_step\":$BEST"
  MODEL=$O/model_$(printf %06d $BEST).pt
  if [ ! -f $O/DEV_DONE ]; then
    $PY $S41/m41_select.py dev_start $O $LABEL; rc=$?
    if [ $rc -eq 11 ]; then ev dev_refused; exit 11; fi
    if [ $rc -eq 0 ]; then
      ev dev_start ",\"best_step\":$BEST"
      D=$O/devM41_$LABEL
      infer $D $MODEL $DEVR || { ev dev_infer_failed; exit 9; }
      ( cd $EXPS/$ENC_EXP/scripts && taskset -c 18-23 $PY c1_eval.py $D M41pc2_${LABEL}_step$BEST DEV --compare-dev >> $D/eval.log 2>&1 && \
        taskset -c 18-23 $PY c1_boot.py $D $D/bootstrap_dev.json >> $D/eval.log 2>&1 ) || { ev dev_eval_failed; exit 9; }
      ( cd $S34 && M34_BACKEND=ntv2 taskset -c 18-23 $GPY $S41/m41_compare.py $O/dev_compare.json $LABEL=$D T2=$O34/pc2_T2/devM34_T2 \
          T2_s2=$O34/pc2_T2_s2/devM34_T2_s2 T2_s3=$O34/pc2_T2_s3/devM34_T2_s3 >> $D/eval.log 2>&1 ) || ev dev_compare_failed
      rm -rf $D/pot
      touch $O/DEV_DONE
      ev dev_done
    fi
  fi
  if [ ! -f $O/setB/SETB_DONE ]; then
    infer $O/setB $MODEL $REGB --setb || { ev setb_infer_failed; exit 9; }
    ( cd $S34 && M34_BACKEND=ntv2 taskset -c 18-23 nice -n 10 $GPY m34_p2b_score.py --new $LABEL=$O/setB --old T2=$O34/p2b/T2 \
        --extra T2_s2=$O34/pc2_T2_s2/setB --extra T2_s3=$O34/pc2_T2_s3/setB --out $O/setB_${LABEL}_result.json --label M41_setB >> $O/setB.log 2>&1 ) || ev setb_score_failed
    rm -rf $O/setB/pot
    echo "{\"arm\":\"$LABEL\",\"step\":$BEST,\"time\":\"$(date -Is)\"}" > $O/setB/SETB_DONE
    ev setb_done
  fi
  touch $O/ARM_DONE
fi
ev arm_done
