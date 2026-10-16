#!/usr/bin/env bash
# CK1 smoke (M6): 1 seed, 2 epochs, --save-raw-pred. Confirm raw_pred npz (val+test) + val_eval_subsets written.
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"; cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh; set +u; conda activate coding-rna; set -u
python -m src.foundation_probe.train_probe_head \
  --cache outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species \
  --species data/m1_screen/saccharomyces_cerevisiae data/m1_screen/drosophila_melanogaster \
  --exp-id FP-FRAGFIX-RP-SMOKE --out-dir outputs/FP-FRAGFIX-RP-SMOKE \
  --head convlstm --window 2048 --sample-fraction 0.3 --epochs 2 --patience 2 --batch-size 16 --lr 1e-3 \
  --seed 0 --class-weighting sqrt_inv --loss fp_aware --fp-lambda 1.0 --postproc constrained --save-raw-pred 2>&1 \
  | grep -ivE 'rematerialization|bfc_allocator|Current allocation|captured constants' | tail -6
python scripts/_verify_m6_rawpred.py
echo "M6_RAWPRED_SMOKE_DONE"
