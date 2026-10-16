#!/usr/bin/env bash
#SBATCH --job-name=M25RDEVREDEC3
#SBATCH --partition=private-teodoro-gpu
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --output=outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R3/logs/%x_%j.out
#SBATCH --error=outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R3/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/home/users/j/jwang/coding-rna"
EXP_ID="M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R3"
SOURCE_ID="M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
STATUS_FILE="${OUT_DIR}/STATUS"
LOCAL_LOG_DIR=""

cd "${ROOT}"
if [[ -d "${OUT_DIR}" ]] && find "${OUT_DIR}" -mindepth 1 -type f ! -path "${OUT_DIR}/logs/*" ! -name JOBID -print -quit | grep -q .; then
  echo "Refusing to overwrite existing Stage 1 output: ${OUT_DIR}" >&2
  exit 3
fi
for epoch in 1 2 3; do
  test -s "${ROOT}/outputs/${SOURCE_ID}/checkpoints/epoch_${epoch}.pt"
done
test -s "${ROOT}/outputs/${SOURCE_ID}/validation_grid_diagnostics.json"

mkdir -p "${OUT_DIR}/logs"
echo "RUNNING" > "${STATUS_FILE}"
copy_python_logs() {
  if [[ -z "${LOCAL_LOG_DIR}" || ! -d "${LOCAL_LOG_DIR}" ]]; then
    return 0
  fi
  local rc=0
  if [[ -f "${LOCAL_LOG_DIR}/python.out" ]]; then
    cp "${LOCAL_LOG_DIR}/python.out" "${OUT_DIR}/logs/python_${SLURM_JOB_ID}.out" || rc=$?
  fi
  if [[ -f "${LOCAL_LOG_DIR}/python.err" ]]; then
    cp "${LOCAL_LOG_DIR}/python.err" "${OUT_DIR}/logs/python_${SLURM_JOB_ID}.err" || rc=$?
  fi
  return "${rc}"
}
finalize_job() {
  local rc=$?
  local copy_rc=0
  trap - EXIT
  copy_python_logs || copy_rc=$?
  if [[ ${rc} -eq 0 && ${copy_rc} -ne 0 ]]; then
    rc=${copy_rc}
  fi
  if [[ ${rc} -ne 0 ]]; then
    echo "FAILED" > "${STATUS_FILE}"
  fi
  exit "${rc}"
}
trap finalize_job EXIT
trap 'exit 143' TERM INT

source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export HF_HOME="${HOME}/.cache/huggingface"
export PYTHONUNBUFFERED=1
LOCAL_LOG_DIR="${SLURM_TMPDIR:-/tmp}/${EXP_ID}_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_LOG_DIR}"

{
  echo "exp_id=${EXP_ID}"
  echo "source_exp_id=${SOURCE_ID}"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "git_commit=$(git rev-parse HEAD)"
  echo "host=$(hostname)"
  date -Is
  python -V
  python - <<'PY'
import torch, transformers, peft
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("transformers", transformers.__version__, "peft", peft.__version__)
PY
} > "${OUT_DIR}/env.txt" 2>&1

set +e
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/redecode_error_decomposition.py \
  --root "${ROOT}" \
  --experiment-id "${EXP_ID}" \
  --config "${ROOT}/configs/${SOURCE_ID}.yaml" \
  --m25r-output "${ROOT}/outputs/${SOURCE_ID}" \
  --out-dir "${OUT_DIR}" \
  > "${LOCAL_LOG_DIR}/python.out" \
  2> "${LOCAL_LOG_DIR}/python.err"
python_rc=$?
set -e
if [[ ${python_rc} -ne 0 ]]; then
  exit "${python_rc}"
fi

test -s "${OUT_DIR}/stage1_diagnostic.json"
test "$(<"${STATUS_FILE}")" = "COMPLETED_STAGE1_REVIEW_REQUIRED"

