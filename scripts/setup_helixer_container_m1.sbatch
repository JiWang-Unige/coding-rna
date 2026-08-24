#!/usr/bin/env bash
#SBATCH --job-name=SETUP-HELIXER-CONTAINER-M1
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/SETUP-HELIXER-CONTAINER-M1/logs/%x_%j.out
#SBATCH --error=outputs/SETUP-HELIXER-CONTAINER-M1/logs/%x_%j.err

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/SETUP-HELIXER-CONTAINER-M1"
LOG_DIR="${OUT_DIR}/logs"
SIF_DIR="${ROOT}/refs/containers"
WEIGHT_DIR="${ROOT}/refs/weights/helixer-2025"
CACHE_DIR="${ROOT}/refs/cache/apptainer"
TMP_DIR="${OUT_DIR}/tmp"
STATUS_FILE="${OUT_DIR}/SETUP_STATUS"
SIF="${SIF_DIR}/helixer-docker_latest.sif"

mkdir -p "${LOG_DIR}" "${SIF_DIR}" "${WEIGHT_DIR}" "${CACHE_DIR}" "${TMP_DIR}"
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

export APPTAINER_CACHEDIR="${CACHE_DIR}"
export APPTAINER_TMPDIR="${TMP_DIR}"
export APPTAINER_BINDPATH="/srv:/srv,/home:/home"

apptainer pull --force "${SIF}" docker://gglyptodon/helixer-docker:latest

apptainer run "${SIF}" fetch_helixer_models.py --lineage fungi --custom-path "${WEIGHT_DIR}"
apptainer run "${SIF}" fetch_helixer_models.py --lineage invertebrate --custom-path "${WEIGHT_DIR}"

find "${SIF_DIR}" "${WEIGHT_DIR}" -maxdepth 3 -type f -printf "%p\t%s bytes\n" | sort > "${OUT_DIR}/artifacts.tsv"

echo "COMPLETED" > "${STATUS_FILE}"
trap - EXIT
