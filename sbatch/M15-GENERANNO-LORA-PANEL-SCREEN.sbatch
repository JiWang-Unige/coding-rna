#!/usr/bin/env bash
#SBATCH --job-name=M15GENPNL
#SBATCH --array=0-1%1
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --time=35:50:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/fp_segnt_logs/%x_%A_%a.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%A_%a.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ARM="${SLURM_ARRAY_TASK_ID:-0}"
SP_LIST=(arabidopsis_thaliana oryza_sativa)

case "${ARM}" in
  0)
    EXP_ID="M15-GENERANNO-1P2B-CDS-PANEL-SCREEN"
    MODEL_NAME="GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
    MODEL_TASK="token_classification"
    WINDOW=6144
    LIMIT_TRAIN=1024
    LIMIT_VAL=512
    ;;
  1)
    EXP_ID="M15-GENERANNO-0P5B-BASE-PANEL-SCREEN"
    MODEL_NAME="GenerTeam/GENERanno-eukaryote-0.5b-base"
    MODEL_TASK="masked_lm"
    WINDOW=1024
    LIMIT_TRAIN=2048
    LIMIT_VAL=1024
    ;;
  *)
    echo "Unknown ARM=${ARM}" >&2
    exit 2
    ;;
esac

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
  echo "model_name=${MODEL_NAME}"
  echo "model_task=${MODEL_TASK}"
  echo "window=${WINDOW}"
  echo "limit_train_windows=${LIMIT_TRAIN}"
  echo "limit_val_windows=${LIMIT_VAL}"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
  echo "host=$(hostname)"
  date -Is
  python -V
  python - <<'PY'
import torch, transformers, peft
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("transformers", transformers.__version__, "peft", peft.__version__)
PY
} > "${OUT_DIR}/env.txt" 2>&1

python -u -m src.foundation_probe.train_generanno_lora_3class \
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --model-name "${MODEL_NAME}" --model-task "${MODEL_TASK}" \
  --window "${WINDOW}" --sample-fraction 0.1 --epochs 2 --patience 1 \
  --batch-size 1 --predict-batch 1 --lr 8e-4 --lora-lr 2e-5 --seed 0 \
  --class-weighting sqrt_inv --loss fp_aware --fp-lambda 1.0 \
  --postproc constrained --min-cds-len 60 --max-fill-gap 20 \
  --lora-r 8 --lora-alpha 16 --lora-dropout 0.05 \
  --lora-targets q_proj,k_proj,v_proj,o_proj --attn-implementation sdpa --bf16 \
  --limit-train-windows "${LIMIT_TRAIN}" --limit-val-windows "${LIMIT_VAL}"

metric_paths=()
for SP in "${SP_LIST[@]}"; do
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${OUT_DIR}/eval_subsets/${SP}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${SP}.gff" \
    --genome-fasta "${OUT_DIR}/eval_subsets/${SP}/genome.fa" \
    --output-json "${METRIC_DIR}/${SP}.metrics.json" \
    --experiment-id "${EXP_ID}_${SP}" --profile screen --span-mode cds
  metric_paths+=("${METRIC_DIR}/${SP}.metrics.json")
done

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${metric_paths[@]}" \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" --profile screen

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
