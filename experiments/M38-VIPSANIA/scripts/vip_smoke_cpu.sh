#!/bin/bash
# M38-VIPSANIA CPU smoke test (GPU hidden: CUDA_VISIBLE_DEVICES="").  A 500-kb slice of the Arabidopsis TRAINING chromosome NC_003070.9 (not DEV, not set B, not panel):
#   1 pretrained annotate (no finetune)            -> smoke/pre/smoke.gff3
#   3 (mode `unmasked`, side check) pretrained annotate of the same slice with the soft-masking stripped (the authors advise stripping it when the masking is of unknown quality)
#   2 finetune plumbing check (1 epoch, 2 micro-steps of 1 window, accumulation 2: tests the code path and the checkpoint, NOT a real finetuning) + annotate -> smoke/ft/smoke.gff3
# Runtime / peak RAM of each step go to smoke/*.res.json (vip_guard.py); a RAM guard kills the group above 40 GB.
# usage: nice -n 19 bash vip_smoke_cpu.sh [pre|ft|both|unmasked]
source "$(dirname "$0")/vip_env.sh"
export CUDA_VISIBLE_DEVICES=""
export TF_FORCE_GPU_ALLOW_GROWTH=true
S=$M38/smoke; mkdir -p $S/pre $S/ft $S/unmasked
CPUS=${VIP_SMOKE_CPUS:-25-31}
GUARD="python3.12 -I $M38_SCRIPTS/vip_guard.py --max-rss-gb ${VIP_SMOKE_MAX_RSS_GB:-40}"
FA=$S/slice.fa
[ -s $FA ] || python3.12 -I $M38_SCRIPTS/vip_slice.py $M38/genomes/arabidopsis_thaliana.nuclear.fa NC_003070.9 5000001 5500000 $FA
what=${1:-both}
if [ "$what" = pre ] || [ "$what" = both ]; then
  $GUARD --out $S/pre.res.json -- nice -n 19 taskset -c $CPUS $VIP annotate $MODEL_ID $FA -o $S/pre/smoke.gff3 -B 1 > $S/pre.out 2> $S/pre.err
  echo "pre rc=$?"
fi
if [ "$what" = ft ] || [ "$what" = both ]; then
  $GUARD --out $S/ft.res.json -- nice -n 19 taskset -c $CPUS $VIP annotate $MODEL_ID $FA -o $S/ft/smoke.gff3 -B 1 --finetune --finetune_B 1 --finetune_epochs 1 \
      -fo trainer.train_steps=2 -fo trainer.gradient_accumulation_steps=2 > $S/ft.out 2> $S/ft.err
  echo "ft rc=$?"
fi
if [ "$what" = unmasked ]; then
  awk '/^>/ {print; next} {print toupper($0)}' $FA > $S/slice_unmasked.fa
  $GUARD --out $S/unmasked.res.json -- nice -n 19 taskset -c $CPUS $VIP annotate $MODEL_ID $S/slice_unmasked.fa -o $S/unmasked/smoke.gff3 -B 1 > $S/unmasked.out 2> $S/unmasked.err
  echo "unmasked rc=$?"
fi
