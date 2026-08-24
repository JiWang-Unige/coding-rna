#!/usr/bin/env bash
#SBATCH --job-name=M8SMOKE
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TB-GBF1-MULTICLASS-M8 CK2 SMOKE — validate the multi-class code path end-to-end on arabidopsis
# (M7-cached features): --label-scheme multiclass (8-class strand-aware) + --decoder crf
# (LinearChainCRFVec 8x8 transitions = biologically-meaningful structure) + fp_aware + constrained.
# Check: runs, per-class(3class-collapsed) macroF1 sane, head not collapsed, predicted_genes>0.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_smoke_M8.sbatch'
set -euo pipefail
EXP_ID="M8-MC-SMOKE"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
CACHE="${ROOT}/outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u

python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP_ARA}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --head convlstm --label-scheme multiclass --decoder crf --window 2048 --sample-fraction 0.2 \
  --epochs 2 --patience 2 --batch-size 16 --lr 1e-3 --seed 0 --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 2>&1 \
  | grep -ivE 'rematerialization|bfc_allocator|Current allocation'

python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/arabidopsis_thaliana.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/genome.fa" \
  --output-json "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" \
  --experiment-id "${EXP_ID}_arabidopsis_thaliana" --profile screen --span-mode cds
python -c "import json;m=json.load(open('${METRIC_DIR}/arabidopsis_thaliana.metrics.json'));print('SMOKE spec',round(m.get('intergenic_specificity',-1),4),'gbF1',round(m.get('gene_body_F1_unconstrained',-1),4),'gcount',round(m.get('predicted_gene_count_ratio_vs_reference',-1),3))"
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "M8_SMOKE_DONE"
