#!/usr/bin/env bash
#SBATCH --job-name=M9SMOKE
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TB-UNFREEZE-BACKBONE-M9 CK2 SANITY SMOKE — validate the unfreeze trainer end-to-end on arabidopsis:
# NT-v2-500m backbone (generanno env) unfreeze top 2 + 3c FP-aware head, tiny budget. Check: loads,
# token alignment, backbone top-layer grad flows, head not collapse (per-class macroF1), no OOM, gbF1 sane.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_smoke_M9.sbatch'
set -euo pipefail
EXP_ID="M9-UNFREEZE-SMOKE"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export HF_HOME="${HOME}/.cache/huggingface"

python -m src.foundation_probe.train_unfreeze_backbone \
  --species "${SP_ARA}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --unfreeze-layers 2 --window 2046 --sample-fraction 0.1 --epochs 2 --patience 2 \
  --batch-size 4 --lr 1e-3 --backbone-lr 1e-5 --seed 0 --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 2>&1 \
  | grep -ivE 'rematerialization|bfc_allocator|Current allocation|Some weights'

python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/arabidopsis_thaliana.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/genome.fa" \
  --output-json "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" \
  --experiment-id "${EXP_ID}_arabidopsis_thaliana" --profile screen --span-mode cds
python -c "import json;m=json.load(open('${METRIC_DIR}/arabidopsis_thaliana.metrics.json'));print('M9 SMOKE spec',round(m.get('intergenic_specificity',-1),4),'gbF1',round(m.get('gene_body_F1_unconstrained',-1),4),'gcount',round(m.get('predicted_gene_count_ratio_vs_reference',-1),3))"
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "M9_SMOKE_DONE"
