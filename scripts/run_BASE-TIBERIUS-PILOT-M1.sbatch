#!/usr/bin/env bash
#SBATCH --job-name=BASE-TIBERIUS-PILOT-M1
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/BASE-TIBERIUS-PILOT-M1/logs/%x_%j.out
#SBATCH --error=outputs/BASE-TIBERIUS-PILOT-M1/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

EXP_ID="BASE-TIBERIUS-PILOT-M1"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
PRED_DIR="${OUT_DIR}/predictions"
METRIC_DIR="${OUT_DIR}/metrics"
LOG_DIR="${OUT_DIR}/logs"
STATUS_FILE="${OUT_DIR}/STATUS"

mkdir -p "${PRED_DIR}" "${METRIC_DIR}" "${LOG_DIR}"
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

export APPTAINER_BINDPATH="/srv:/srv,/home:/home"
export SINGULARITY_BINDPATH="${APPTAINER_BINDPATH}"
export TMPDIR="${SLURM_TMPDIR:-${OUT_DIR}/tmp}"
mkdir -p "${TMPDIR}"

run_tiberius() {
  local slug="$1"
  local model_cfg="$2"
  local genome="$3"
  local reference="$4"
  local pred="${PRED_DIR}/${slug}.gtf"
  local metric="${METRIC_DIR}/${slug}.metrics.json"

  echo "[${EXP_ID}] Running Tiberius ${model_cfg} for ${slug}"
  python refs/repos/tiberius-2024/tiberius.py \
    --singularity \
    --genome "${genome}" \
    --model_cfg "${model_cfg}" \
    --out "${pred}" \
    --batch_size 4 \
    --id_prefix "${slug}_"

  echo "[${EXP_ID}] Evaluating ${slug}"
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${reference}" \
    --prediction-gtf "${pred}" \
    --genome-fasta "${genome}" \
    --output-json "${metric}" \
    --experiment-id "${EXP_ID}_${slug}" \
    --profile screen
}

run_tiberius \
  "saccharomyces_cerevisiae" \
  "fungi" \
  "data/m1_screen/saccharomyces_cerevisiae/genome.fa" \
  "data/m1_screen/saccharomyces_cerevisiae/reference.gff3"

run_tiberius \
  "drosophila_melanogaster" \
  "insecta" \
  "data/m1_screen/drosophila_melanogaster/genome.fa" \
  "data/m1_screen/drosophila_melanogaster/reference.gff3"

python scripts/aggregate_gene_body_metrics.py \
  --metrics \
    "${METRIC_DIR}/saccharomyces_cerevisiae.metrics.json" \
    "${METRIC_DIR}/drosophila_melanogaster.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" \
  --profile screen

set +e
python scripts/validate_goal.py \
  --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" \
  --profile screen \
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
