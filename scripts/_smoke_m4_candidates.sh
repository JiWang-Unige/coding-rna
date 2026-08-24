#!/usr/bin/env bash
# CK1 sanity smoke for TA-FOUNDATION-DECODER-M4: each of the 3 candidates, 1 seed, 2 epochs,
# on the FULL FEATCACHE. Confirms the new code paths (fp_aware loss / raw-DNA fusion / CRF) run
# end-to-end and DON'T collapse (per_class non-degenerate), before the 15-run batch.
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u
CACHE="outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
SP_SAC="data/m1_screen/saccharomyces_cerevisiae"
SP_DMEL="data/m1_screen/drosophila_melanogaster"

run_one() {
  local CAND="$1"; shift
  local EXP="FP-SEGNT-${CAND}-SMOKE"; local OUT="outputs/${EXP}"; local MET="${OUT}/metrics"
  mkdir -p "${MET}"
  echo "########## SMOKE ${CAND} ($*) ##########"
  python -m src.foundation_probe.train_probe_head \
    --cache "${CACHE}" --species "${SP_SAC}" "${SP_DMEL}" \
    --exp-id "${EXP}" --out-dir "${OUT}" \
    --head convlstm --window 2048 --sample-fraction 0.3 --epochs 2 --patience 2 \
    --batch-size 16 --lr 1e-3 --seed 0 --class-weighting sqrt_inv "$@" 2>&1 \
    | grep -ivE 'rematerialization|bfc_allocator|Current allocation|^W[0-9]|captured constants'
  for sp in saccharomyces_cerevisiae drosophila_melanogaster; do
    python scripts/eval_gene_body_mask.py \
      --reference-gtf "${OUT}/eval_subsets/${sp}/reference.gff3" \
      --prediction-gtf "${OUT}/predictions/${sp}.gff" \
      --genome-fasta "${OUT}/eval_subsets/${sp}/genome.fa" \
      --output-json "${MET}/${sp}.metrics.json" \
      --experiment-id "${EXP}_${sp}" --profile screen --span-mode cds >/dev/null 2>&1
  done
  python scripts/aggregate_gene_body_metrics.py \
    --metrics "${MET}/saccharomyces_cerevisiae.metrics.json" "${MET}/drosophila_melanogaster.metrics.json" \
    --output-json "${MET}/metrics.json" --experiment-id "${EXP}" --profile screen >/dev/null 2>&1
  python -c "import json;d=json.load(open('${MET}/metrics.json'));print('  ${CAND} SMOKE:',{k:(round(d[k],4) if isinstance(d.get(k),float) else d.get(k)) for k in ('intergenic_specificity','macro_intergenic_specificity','gene_body_F1_unconstrained','intergenic_FPR','predicted_gene_count_ratio_vs_reference')})"
}

run_one FPLOSS --loss fp_aware --fp-lambda 1.0
run_one FUSION --fuse-raw-dna
run_one CRF    --decoder crf --loss fp_aware --fp-lambda 1.0 --crf-aux-ce 1.0
echo "ALL_M4_SMOKES_DONE"
