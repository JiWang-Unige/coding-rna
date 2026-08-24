#!/usr/bin/env bash
#SBATCH --job-name=M9ARM
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=23:50:00
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TB-UNFREEZE-BACKBONE-M9 CK3 (single-species fix). arabidopsis only, epochs=4.
# N=0 frozen CONTROL; N=2/4 staged unfreeze. generanno env. NON-CLAIM.
set -euo pipefail
N="${1:?usage: <unfreeze 0|2|4> <seed>}"; SEED="${2:?seed}"
EXP_ID="M9-UNFREEZE-L${N}-s${SEED}"
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
export PYTHONUNBUFFERED=1

stdbuf -oL -eL python -u -m src.foundation_probe.train_unfreeze_backbone \
  --species "${SP_ARA}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --unfreeze-layers "${N}" --window 2046 --sample-fraction 0.3 --epochs 4 --patience 2 \
  --batch-size 4 --lr 1e-3 --backbone-lr 1e-5 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 2>&1 \
  | stdbuf -oL grep -ivE "rematerialization|bfc_allocator|Current allocation|Some weights"

python scripts/eval_gene_body_mask.py --reference-gtf "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/arabidopsis_thaliana.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/arabidopsis_thaliana/genome.fa" \
  --output-json "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" \
  --experiment-id "${EXP_ID}_arabidopsis_thaliana" --profile screen --span-mode cds
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
