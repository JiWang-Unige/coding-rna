#!/usr/bin/env bash
#SBATCH --job-name=M16MULTI
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --time=71:50:00
#SBATCH --partition=shared-gpu
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
SEED=0
EXP_ID="M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s${SEED}"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
STATUS_FILE="${OUT_DIR}/STATUS"

mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
python "${ROOT}/scripts/pre_submit_gate.py" --exp-id "${EXP_ID}" "${ROOT}" --format text

if [[ -e "${STATUS_FILE}" ]] && grep -qx "COMPLETED" "${STATUS_FILE}"; then
  echo "${EXP_ID} already COMPLETED; refusing to overwrite."
  exit 0
fi
if [[ -d "${OUT_DIR}" ]] && find "${OUT_DIR}" -mindepth 1 -type f ! -name code_review_gate.json ! -name STATUS -print -quit | grep -q .; then
  echo "Refusing to overwrite non-empty output dir: ${OUT_DIR}" >&2
  exit 3
fi

echo "RUNNING" > "${STATUS_FILE}"
mark_failed() {
  rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "FAILED" > "${STATUS_FILE}"
  fi
  exit ${rc}
}
trap mark_failed EXIT

cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export HF_HOME="${HOME}/.cache/huggingface"
export PYTHONUNBUFFERED=1

{
  echo "exp_id=${EXP_ID}"
  echo "protocol=Arabidopsis+rice train+val calibration -> Arabidopsis lyrata + gallus + drosophila tests"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "host=$(hostname)"
  date -Is
  python -V
  python - <<'PY'
import torch, transformers
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("transformers", transformers.__version__)
PY
} > "${OUT_DIR}/env.txt" 2>&1

python -u -m src.foundation_probe.train_unfreeze_backbone \
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa data/m1_screen/arabidopsis_lyrata data/m1_screen/gallus_gallus data/m1_screen/drosophila_melanogaster \
  --train-species-names arabidopsis_thaliana oryza_sativa \
  --val-species-names arabidopsis_thaliana oryza_sativa \
  --test-species-names arabidopsis_lyrata gallus_gallus drosophila_melanogaster \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --unfreeze-layers 12 --window 2046 --sample-fraction 0.3 --epochs 4 --patience 2 \
  --batch-size 4 --lr 1e-3 --backbone-lr 1e-5 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 \
  --predict-batch 8 --save-raw-scores

python -u scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --species arabidopsis_thaliana oryza_sativa arabidopsis_lyrata gallus_gallus drosophila_melanogaster \
  --val-species arabidopsis_thaliana oryza_sativa \
  --test-species arabidopsis_lyrata gallus_gallus drosophila_melanogaster \
  --profile screen \
  --intergenic-biases 0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0 \
  --min-cds-lens 60,90,120 --max-fill-gaps 0,20 \
  --target-fpr 0.01 --min-gbf1 0.50 --max-gene-count-ratio 1.25

echo "COMPLETED" > "${STATUS_FILE}"

set +e
python scripts/validate_goal.py \
  --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" \
  --profile screen \
  --run-status "${STATUS_FILE}" > "${METRIC_DIR}/validate_goal.json"
VALIDATE_RC=$?
set -e
echo "${VALIDATE_RC}" > "${METRIC_DIR}/validate_goal.exit_code"
if [[ ${VALIDATE_RC} -eq 3 ]]; then
  echo "FAILED" > "${STATUS_FILE}"
  exit 3
fi

trap - EXIT
echo "DONE ${EXP_ID}"
