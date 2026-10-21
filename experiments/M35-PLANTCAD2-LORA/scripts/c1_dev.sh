#!/bin/bash
# one-shot DEV evaluation of the selected checkpoint (run exactly once, after c1_stage.sh): inference on both DEV chromosomes + same-ruler evaluation with comparison rows.
# usage: c1_dev.sh <runname>
RUN=$1
O=$HOME/coding-rna/outputs/M35-PLANTCAD2-LORA; R=$O/$RUN
BEST=$(cat $R/best_step.txt)
PY=$HOME/envs/plantcad2/bin/python
cd $HOME/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA/scripts
if [ -e $R/dev_$BEST/DEV_STARTED ]; then echo "DEV already started for $BEST: refusing a second DEV run"; exit 3; fi
mkdir -p $R/dev_$BEST && date -Is > $R/dev_$BEST/DEV_STARTED
bash c1_infer_retry.sh $R/dev_$BEST $R/model_$(printf %06d $BEST).pt NC_003074.8,NC_089041.1 4
taskset -c ${M35_SCORE_CPUS:-18-23} $PY c1_eval.py $R/dev_$BEST step$BEST DEV --compare-dev | tee $R/dev_eval_step$BEST.log
rm -rf $R/dev_$BEST/pot
