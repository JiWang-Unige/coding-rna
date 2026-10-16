#!/usr/bin/env bash
set -euo pipefail

# Run validation-only decode/FPR calibration for completed M19 seed outputs.
# This intentionally selects decode parameters on VAL only, then applies them once to TEST.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$ROOT"

species=(arabidopsis_thaliana oryza_sativa)

for seed in 0 1; do
  exp_id="M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s${seed}"
  out_dir="outputs/${exp_id}"
  missing=0
  for sp in "${species[@]}"; do
    for split in val test; do
      if [[ ! -s "${out_dir}/raw_scores/${split}_${sp}.npz" ]]; then
        echo "[skip] ${exp_id}: missing ${out_dir}/raw_scores/${split}_${sp}.npz" >&2
        missing=1
      fi
    done
  done
  if [[ "$missing" -eq 1 ]]; then
    continue
  fi

  # calibrate_decode.py writes final TEST predictions/metrics back to the run root.
  # Keep the original raw-decode artifacts so M19 can report both operating points.
  backup_dir="${out_dir}/metrics/pre_rawcal"
  if [[ -s "${out_dir}/metrics/metrics.json" && ! -s "${backup_dir}/metrics.json" ]]; then
    mkdir -p "${backup_dir}"
    cp -a "${out_dir}/metrics/"*.json "${backup_dir}/"
  fi
  if [[ -d "${out_dir}/predictions" && ! -d "${out_dir}/predictions_pre_rawcal" ]]; then
    cp -a "${out_dir}/predictions" "${out_dir}/predictions_pre_rawcal"
  fi

  python3 scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py \
    --exp-id "${exp_id}" \
    --out-dir "${out_dir}" \
    --species "${species[@]}" \
    --profile screen \
    --intergenic-biases "0,0.25,0.5,0.75,1.0,1.25,1.5,2.0" \
    --min-cds-lens "60,90,120" \
    --max-fill-gaps "0,20" \
    --target-fpr 0.01 \
    --min-gbf1 0.70 \
    --max-gene-count-ratio 1.25
done
