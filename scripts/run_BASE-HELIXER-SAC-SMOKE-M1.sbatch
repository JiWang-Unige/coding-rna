#!/usr/bin/env bash
#SBATCH --job-name=BASE-HELIXER-SAC-SMOKE-M1
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/BASE-HELIXER-SAC-SMOKE-M1/logs/%x_%j.out
#SBATCH --error=outputs/BASE-HELIXER-SAC-SMOKE-M1/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

EXP_ID="BASE-HELIXER-SAC-SMOKE-M1"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
PRED_DIR="${OUT_DIR}/predictions"
METRIC_DIR="${OUT_DIR}/metrics"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
STATUS_FILE="${OUT_DIR}/STATUS"
SIF="${ROOT}/refs/containers/helixer-docker_latest.sif"
WEIGHT_DIR="${ROOT}/refs/weights/helixer-2025"

mkdir -p "${PRED_DIR}" "${METRIC_DIR}" "${LOG_DIR}" "${TMP_DIR}"
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
conda activate coding-rna

if [[ ! -s "${SIF}" ]]; then
  echo "Missing Helixer SIF: ${SIF}" >&2
  exit 2
fi
if [[ ! -s "${WEIGHT_DIR}/fungi/fungi_v0.3_a_0100.h5" ]]; then
  echo "Missing Helixer fungi model under ${WEIGHT_DIR}" >&2
  exit 2
fi

export APPTAINER_BINDPATH="/srv:/srv,/home:/home"
export APPTAINERENV_CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"

apptainer run --nv "${SIF}" Helixer.py \
  --lineage fungi \
  --downloaded-model-path "${WEIGHT_DIR}" \
  --fasta-path "${ROOT}/data/m1_screen/saccharomyces_cerevisiae/genome.fa" \
  --species "Saccharomyces_cerevisiae" \
  --temporary-dir "${TMP_DIR}" \
  --batch-size 16 \
  --gff-output-path "${PRED_DIR}/saccharomyces_cerevisiae.gff3"

python scripts/eval_gene_body_mask.py \
  --reference-gtf data/m1_screen/saccharomyces_cerevisiae/reference.gff3 \
  --prediction-gtf "${PRED_DIR}/saccharomyces_cerevisiae.gff3" \
  --genome-fasta data/m1_screen/saccharomyces_cerevisiae/genome.fa \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" \
  --profile smoke

set +e
python scripts/validate_goal.py \
  --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" \
  --profile smoke \
  --run-status "${STATUS_FILE}" \
  > "${METRIC_DIR}/validate_goal.json"
VALIDATE_RC=$?
set -e
echo "${VALIDATE_RC}" > "${METRIC_DIR}/validate_goal.exit_code"
if [[ ${VALIDATE_RC} -eq 3 ]]; then
  echo "FAILED" > "${STATUS_FILE}"
  exit 3
fi

echo "COMPLETED" > "${STATUS_FILE}"
trap - EXIT
