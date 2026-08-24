#!/usr/bin/env bash
#SBATCH --job-name=M12BEXT
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_a100-pcie-40gb:1
#SBATCH --time=11:50:00
#SBATCH --partition=shared-gpu
#SBATCH --output=outputs/M12B-SAMEPANEL-BASELINES/logs/%x_%A_%a.out
#SBATCH --error=outputs/M12B-SAMEPANEL-BASELINES/logs/%x_%A_%a.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ARMS=(ANNEVO TIBERIUS HELIXER)
ARM="${ARMS[${SLURM_ARRAY_TASK_ID:-0}]}"
EXP_ID="M12B-SAMEPANEL-BASELINES-${ARM}"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
PRED_DIR="${OUT_DIR}/predictions"
METRIC_DIR="${OUT_DIR}/metrics"
LOG_DIR="${ROOT}/outputs/M12B-SAMEPANEL-BASELINES/logs"
TMP_DIR="/tmp/m12b_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
STATUS_FILE="${OUT_DIR}/STATUS"
SPECIES=(arabidopsis_thaliana oryza_sativa)

mkdir -p "${PRED_DIR}" "${METRIC_DIR}" "${LOG_DIR}" "${TMP_DIR}"
python "${ROOT}/scripts/pre_submit_gate.py" --exp-id "${EXP_ID}" "${ROOT}" --format text

if [[ -e "${STATUS_FILE}" ]] && grep -qx "COMPLETED" "${STATUS_FILE}"; then
  echo "${EXP_ID} already COMPLETED; refusing to overwrite."
  exit 0
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
export PYTHONUNBUFFERED=1
export APPTAINER_BINDPATH="/srv:/srv,/home:/home"
export SINGULARITY_BINDPATH="${APPTAINER_BINDPATH}"
export APPTAINERENV_CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export TMPDIR="${TMP_DIR}/jobtmp"
mkdir -p "${TMPDIR}"

{
  echo "exp_id=${EXP_ID}"
  echo "arm=${ARM}"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
  echo "host=$(hostname)"
  date -Is
} > "${OUT_DIR}/env.txt"

run_eval_one() {
  local sp="$1"
  local pred="$2"
  local metric="${METRIC_DIR}/${sp}.metrics.json"
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "data/m1_screen/${sp}/reference.gff3" \
    --prediction-gtf "${pred}" \
    --genome-fasta "data/m1_screen/${sp}/genome.fa" \
    --output-json "${metric}" \
    --experiment-id "${EXP_ID}_${sp}" \
    --profile screen --span-mode cds
}

run_annevo() {
  set +u
  conda activate annevo
  set -u
  python - <<'PY' >> "${OUT_DIR}/env.txt" 2>&1
import torch
print("python_arm=annevo")
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
PY
  local annevo_dir="${ROOT}/refs/repos/annevo-2026"
  for sp in "${SPECIES[@]}"; do
    local pred="${PRED_DIR}/${sp}.gff"
    echo "[${EXP_ID}] ANNEVO Magnoliopsida ${sp}"
    ( cd "${annevo_dir}" && python annotation.py \
        -g "${ROOT}/data/m1_screen/${sp}/genome.fa" \
        -m "saved_model/ANNEVO_Magnoliopsida.pt" \
        -l "Magnoliopsida" \
        -o "${pred}" \
        --batch_size 32 -t "${SLURM_CPUS_PER_TASK:-16}" --show_log )
    run_eval_one "${sp}" "${pred}"
  done
}

run_tiberius() {
  set +u
  conda activate coding-rna
  set -u
  python - <<'PY' >> "${OUT_DIR}/env.txt" 2>&1
import sys
print("python_arm=tiberius")
print(sys.version)
PY
  for sp in "${SPECIES[@]}"; do
    local pred="${PRED_DIR}/${sp}.gtf"
    echo "[${EXP_ID}] Tiberius angiosperms ${sp}"
    python refs/repos/tiberius-2024/tiberius.py \
      --singularity \
      --genome "data/m1_screen/${sp}/genome.fa" \
      --model_cfg "angiosperms" \
      --out "${pred}" \
      --batch_size 4 \
      --id_prefix "${sp}_"
    run_eval_one "${sp}" "${pred}"
  done
}

run_helixer() {
  set +u
  conda activate coding-rna
  set -u
  local sif="${ROOT}/refs/containers/helixer-docker_latest.sif"
  local model="${ROOT}/refs/weights/helixer-2025/land_plant/land_plant_v0.3_m_0100.h5"
  [[ -s "${sif}" ]] || { echo "Missing Helixer SIF: ${sif}" >&2; exit 2; }
  [[ -s "${model}" ]] || { echo "Missing Helixer model: ${model}" >&2; exit 2; }
  for sp in "${SPECIES[@]}"; do
    local pred="${PRED_DIR}/${sp}.gff3"
    local sp_tmp="${TMP_DIR}/${sp}"
    mkdir -p "${sp_tmp}"
    echo "[${EXP_ID}] Helixer land_plant ${sp}"
    apptainer run --nv "${sif}" Helixer.py \
      --lineage land_plant \
      --model-filepath "${model}" \
      --subsequence-length 64152 \
      --fasta-path "${ROOT}/data/m1_screen/${sp}/genome.fa" \
      --species "${sp}" \
      --temporary-dir "${sp_tmp}" \
      --batch-size 16 \
      --gff-output-path "${pred}"
    run_eval_one "${sp}" "${pred}"
  done
}

case "${ARM}" in
  ANNEVO) run_annevo ;;
  TIBERIUS) run_tiberius ;;
  HELIXER) run_helixer ;;
  *) echo "Unknown ARM=${ARM}" >&2; exit 2 ;;
esac

metric_args=()
for sp in "${SPECIES[@]}"; do
  metric_args+=("${METRIC_DIR}/${sp}.metrics.json")
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${metric_args[@]}" \
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
