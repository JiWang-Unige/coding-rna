#!/usr/bin/env bash
# CK3 smoke: validate the full probe chain on the YEAST-ONLY cache (exists while fly extracts).
# train light head (2 epochs) -> predict -> CDS GFF -> eval --span-mode cds -> aggregate -> validate.
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u
CACHE="outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
EXP="FP-SEGNT-SMOKE-yeast-convlstm"
OUT="outputs/${EXP}"; MET="${OUT}/metrics"
mkdir -p "${MET}"
SP="data/m1_screen/saccharomyces_cerevisiae"

python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP}" \
  --exp-id "${EXP}" --out-dir "${OUT}" \
  --head convlstm --window 2048 --sample-fraction 0.3 --epochs 3 --patience 2 \
  --batch-size 16 --lr 1e-3 --seed 0 --class-weighting sqrt_inv

python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT}/eval_subsets/saccharomyces_cerevisiae/reference.gff3" \
  --prediction-gtf "${OUT}/predictions/saccharomyces_cerevisiae.gff" \
  --genome-fasta "${OUT}/eval_subsets/saccharomyces_cerevisiae/genome.fa" \
  --output-json "${MET}/saccharomyces_cerevisiae.metrics.json" \
  --experiment-id "${EXP}_sac" --profile screen --span-mode cds

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${MET}/saccharomyces_cerevisiae.metrics.json" \
  --output-json "${MET}/metrics.json" --experiment-id "${EXP}" --profile screen

echo "=== SMOKE metrics (yeast-only) ==="
python -c "import json;d=json.load(open('${MET}/metrics.json'));print({k:d.get(k) for k in ('primary_metric','intergenic_specificity','macro_intergenic_specificity','gene_body_F1_unconstrained','intergenic_FPR','predicted_gene_count_ratio_vs_reference')})"
echo "SMOKE_CHAIN_OK"
