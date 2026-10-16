#!/usr/bin/env bash
# FP-SEGMENTNT-PROBE-M1 feature extraction for both pilot species into a shared cache.
# Run via srun on a GPU node. Usage: bash scripts/_extract_fp_featcache.sh
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate coding-rna
set -u
export HF_HOME="${HOME}/.cache/huggingface"
CACHE="outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
mkdir -p "${CACHE}"
for SP in saccharomyces_cerevisiae drosophila_melanogaster; do
  echo ">>> extracting ${SP}"
  python -m src.foundation_probe.extract_segmentnt \
    --species "data/m1_screen/${SP}" --out-dir "${CACHE}" \
    --model segment_nt_multi_species --tile-tokens 1000 2>&1 \
    | grep -ivE 'rematerialization|bfc_allocator|Current allocation|^W[0-9]|captured constants' \
    | tail -25
done
echo "ALL_EXTRACT_DONE"
