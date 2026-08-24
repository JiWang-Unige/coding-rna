#!/usr/bin/env bash
#SBATCH --job-name=M11L12CAL
#SBATCH --array=0-2
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --time=71:50:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/fp_segnt_logs/%x_%A_%a.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%A_%a.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
SEEDS=(0 1 2)
SEED="${SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}"
EXP_ID="M11-L12-SPEC-CALIBRATION-s${SEED}"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
STATUS_FILE="${OUT_DIR}/STATUS"
SPECIES=(arabidopsis_thaliana oryza_sativa)

mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() {
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "FAILED" > "${STATUS_FILE}"
  fi
  exit $rc
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
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
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
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --unfreeze-layers 12 --window 2046 --sample-fraction 0.3 --epochs 4 --patience 2 \
  --batch-size 4 --lr 1e-3 --backbone-lr 1e-5 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 \
  --predict-batch 8 --save-raw-scores

python -u scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" --species "${SPECIES[@]}" \
  --profile screen \
  --intergenic-biases 0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0 \
  --min-cds-lens 60,90,120 --max-fill-gaps 0,20 \
  --target-fpr 0.01 --min-gbf1 0.70 --max-gene-count-ratio 1.25

echo "COMPLETED" > "${STATUS_FILE}"

python scripts/validate_goal.py \
  --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" \
  --profile screen \
  --run-status "${STATUS_FILE}" > "${METRIC_DIR}/validate_goal.json" || true

trap - EXIT
echo "DONE ${EXP_ID}"
