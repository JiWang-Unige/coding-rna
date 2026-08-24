#!/usr/bin/env bash
#SBATCH --job-name=BASE-ANNEVO-SAC-DMEL-SMOKE-M1
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/BASE-ANNEVO-SAC-DMEL-SMOKE-M1/logs/%x_%j.out
#SBATCH --error=outputs/BASE-ANNEVO-SAC-DMEL-SMOKE-M1/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

# BASE-ANNEVO-SAC-DMEL-SMOKE-M1: ANNEVO ("ANNEVO-light") two-species smoke — the third
# full gene-caller baseline (after Tiberius + Helixer) toward screen_anchor.
# Runs in the dedicated `annevo` conda env (Py3.10/torch2.1/cu12.1; sanctioned exception to
# the coding-rna project env — ANNEVO deps are incompatible). Eval/aggregate/validate use
# only the stdlib so they run fine under py3.10. Eval uses --span-mode cds (the fair
# cross-tool gene-body span per docs/11). Profile=screen; cannot claim.

EXP_ID="BASE-ANNEVO-SAC-DMEL-SMOKE-M1"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ANNEVO_DIR="${ROOT}/refs/repos/annevo-2026"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
PRED_DIR="${OUT_DIR}/predictions"
METRIC_DIR="${OUT_DIR}/metrics"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
STATUS_FILE="${OUT_DIR}/STATUS"

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

source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u                      # annevo env's conda activate.d (libblas_mkl) references unbound MKL_* vars
conda activate annevo
set -u
# ANNEVO decoding uses multiprocessing (-t); its pymp Unix-socket path = $TMPDIR/pymp-*/...
# must stay under the ~108-char AF_UNIX limit. The beegfs output tmp path is too long
# (-> "OSError: AF_UNIX path too long"), so use a SHORT node-local /tmp dir.
export TMPDIR="/tmp/annevo_${SLURM_JOB_ID:-$$}"
mkdir -p "${TMPDIR}"

python -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda, 'avail', torch.cuda.is_available())"

run_annevo() {
  local slug="$1"
  local clade="$2"      # ANNEVO model clade (Fungi / Insecta / ...)
  local lineage="$3"    # ANNEVO -l lineage
  local genome="$4"
  local reference="$5"
  local pred="${PRED_DIR}/${slug}.gff"
  local metric="${METRIC_DIR}/${slug}.metrics.json"

  echo "[${EXP_ID}] Running ANNEVO clade=${clade} lineage=${lineage} for ${slug}"
  ( cd "${ANNEVO_DIR}" && python annotation.py \
      -g "${ROOT}/${genome}" \
      -m "saved_model/ANNEVO_${clade}.pt" \
      -l "${lineage}" \
      -o "${pred}" \
      --batch_size 32 \
      -t "${SLURM_CPUS_PER_TASK:-16}" \
      --show_log )

  echo "[${EXP_ID}] Evaluating ${slug} (--span-mode cds)"
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${reference}" \
    --prediction-gtf "${pred}" \
    --genome-fasta "${genome}" \
    --output-json "${metric}" \
    --experiment-id "${EXP_ID}_${slug}" \
    --profile screen --span-mode cds
}

cd "${ROOT}"

run_annevo \
  "saccharomyces_cerevisiae" \
  "Fungi" "Fungi" \
  "data/m1_screen/saccharomyces_cerevisiae/genome.fa" \
  "data/m1_screen/saccharomyces_cerevisiae/reference.gff3"

run_annevo \
  "drosophila_melanogaster" \
  "Insecta" "Insecta" \
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
