#!/usr/bin/env bash
# CK1 smoke (M5): FPLOSS + --postproc constrained, 1 seed, 2 epochs. Confirm it runs + cuts
# gene_count_ratio (vs FPLOSS 2.25) without collapse.
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"; cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh; set +u; conda activate coding-rna; set -u
CACHE="outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
SP_SAC="data/m1_screen/saccharomyces_cerevisiae"; SP_DMEL="data/m1_screen/drosophila_melanogaster"
EXP="FP-FRAGFIX-CONSTR-SMOKE"; OUT="outputs/${EXP}"; MET="${OUT}/metrics"; mkdir -p "${MET}"
python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP_SAC}" "${SP_DMEL}" --exp-id "${EXP}" --out-dir "${OUT}" \
  --head convlstm --window 2048 --sample-fraction 0.3 --epochs 2 --patience 2 \
  --batch-size 16 --lr 1e-3 --seed 0 --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained 2>&1 \
  | grep -ivE 'rematerialization|bfc_allocator|Current allocation|^W[0-9]|captured constants'
for sp in saccharomyces_cerevisiae drosophila_melanogaster; do
  python scripts/eval_gene_body_mask.py --reference-gtf "${OUT}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT}/predictions/${sp}.gff" --genome-fasta "${OUT}/eval_subsets/${sp}/genome.fa" \
    --output-json "${MET}/${sp}.metrics.json" --experiment-id "${EXP}_${sp}" --profile screen --span-mode cds >/dev/null 2>&1
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${MET}/saccharomyces_cerevisiae.metrics.json" "${MET}/drosophila_melanogaster.metrics.json" \
  --output-json "${MET}/metrics.json" --experiment-id "${EXP}" --profile screen >/dev/null 2>&1
python -c "import json;d=json.load(open('${MET}/metrics.json'));print('CONSTR SMOKE:',{k:(round(d[k],4) if isinstance(d.get(k),float) else d.get(k)) for k in ('intergenic_specificity','macro_intergenic_specificity','gene_body_F1_unconstrained','intergenic_FPR','predicted_gene_count_ratio_vs_reference')})"
echo "(FPLOSS-no-postproc baseline was: spec 0.983 / gbF1 0.578 / gene_count_ratio 2.13 at 2ep)"
echo "M5_CONSTR_SMOKE_DONE"
