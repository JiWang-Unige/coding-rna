#!/usr/bin/env bash
#SBATCH --job-name=BASE-HELIXER-SAC-DMEL-SMOKE-M1
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/logs/%x_%j.out
#SBATCH --error=outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

# BASE-HELIXER-SAC-DMEL-SMOKE-M1: Helixer two-species smoke under the frozen M1 evaluator.
# Mirrors the proven BASE-TIBERIUS-PILOT-M1 two-species + aggregate pattern, swapping the
# caller for Helixer (fungi model for S. cerevisiae, invertebrate model for D. melanogaster).
# Profile=smoke: validates the Helixer runner/evaluator path; cannot claim, cannot set
# screen_anchor. Reports base-weighted + macro + per-species (aggregate_gene_body_metrics.py).

EXP_ID="BASE-HELIXER-SAC-DMEL-SMOKE-M1"
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
  echo "Missing Helixer fungi model under ${WEIGHT_DIR}/fungi" >&2
  exit 2
fi
if [[ ! -s "${WEIGHT_DIR}/invertebrate/invertebrate_v0.3_m_0100.h5" ]]; then
  echo "Missing Helixer invertebrate model under ${WEIGHT_DIR}/invertebrate" >&2
  exit 2
fi

export APPTAINER_BINDPATH="/srv:/srv,/home:/home"
export SINGULARITY_BINDPATH="${APPTAINER_BINDPATH}"
export APPTAINERENV_CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"

run_helixer() {
  local slug="$1"
  local lineage="$2"
  local model_file="$3"
  local subseq_len="$4"
  local species="$5"
  local genome="$6"
  local reference="$7"
  local pred="${PRED_DIR}/${slug}.gff3"
  local metric="${METRIC_DIR}/${slug}.metrics.json"
  local sp_tmp="${TMP_DIR}/${slug}"
  mkdir -p "${sp_tmp}"

  # --model-filepath pins the EXACT downloaded .h5, bypassing the default lineage
  # resolution and its "check for newest model" network call — robust on offline
  # compute nodes (private-teodoro-gpu nodes may have no internet). Manual model
  # selection REQUIRES --subsequence-length (Helixer recommends 21384 fungi /
  # 64152 plants / 213840 animals); we pass the lineage-appropriate value.
  echo "[${EXP_ID}] Running Helixer lineage=${lineage} model=${model_file} subseq=${subseq_len} for ${slug}"
  apptainer run --nv "${SIF}" Helixer.py \
    --lineage "${lineage}" \
    --model-filepath "${WEIGHT_DIR}/${model_file}" \
    --subsequence-length "${subseq_len}" \
    --fasta-path "${ROOT}/${genome}" \
    --species "${species}" \
    --temporary-dir "${sp_tmp}" \
    --batch-size 16 \
    --gff-output-path "${pred}"

  echo "[${EXP_ID}] Evaluating ${slug}"
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${reference}" \
    --prediction-gtf "${pred}" \
    --genome-fasta "${genome}" \
    --output-json "${metric}" \
    --experiment-id "${EXP_ID}_${slug}" \
    --profile smoke
}

run_helixer \
  "saccharomyces_cerevisiae" \
  "fungi" \
  "fungi/fungi_v0.3_a_0100.h5" \
  "21384" \
  "Saccharomyces_cerevisiae" \
  "data/m1_screen/saccharomyces_cerevisiae/genome.fa" \
  "data/m1_screen/saccharomyces_cerevisiae/reference.gff3"

run_helixer \
  "drosophila_melanogaster" \
  "invertebrate" \
  "invertebrate/invertebrate_v0.3_m_0100.h5" \
  "213840" \
  "Drosophila_melanogaster" \
  "data/m1_screen/drosophila_melanogaster/genome.fa" \
  "data/m1_screen/drosophila_melanogaster/reference.gff3"

python scripts/aggregate_gene_body_metrics.py \
  --metrics \
    "${METRIC_DIR}/saccharomyces_cerevisiae.metrics.json" \
    "${METRIC_DIR}/drosophila_melanogaster.metrics.json" \
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
# smoke runner validation: a completed_poor / not_yet baseline is a VALID smoke outcome.
# Only a true failed_run (exit 3: degenerate/infra) marks STATUS=FAILED.
if [[ ${VALIDATE_RC} -eq 3 ]]; then
  echo "FAILED" > "${STATUS_FILE}"
  exit 3
fi

echo "COMPLETED" > "${STATUS_FILE}"
trap - EXIT
