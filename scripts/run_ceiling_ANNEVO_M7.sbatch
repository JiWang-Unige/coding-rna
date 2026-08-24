#!/usr/bin/env bash
#SBATCH --job-name=CEILM7
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=outputs/REANCHOR-CEILING-ANNEVO-M7/logs/%x_%j.out
#SBATCH --error=outputs/REANCHOR-CEILING-ANNEVO-M7/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# REANCHOR-HELDOUT-M7 CK4 — held-out pretrained CEILING via ANNEVO (clade-matched: arabidopsis=
# Magnoliopsida, chicken=Aves). CONTEXTUAL (not the Pareto gate). Eval is restricted to the SAME
# TEST chromosomes the same-budget anchor used (eval_subsets from SCREENREF-tiberius_like-ho-s0)
# + the NEW full-transcript intergenic ruler (--span-mode cds) -> comparable to anchor/candidate.
# Helixer/Tiberius ceiling DEFERRED (weights need download; ANNEVO=published-SOTA candidate is the
# key reference). NON-CLAIM. Usage: ssh baobab 'cd .. && sbatch scripts/run_ceiling_ANNEVO_M7.sbatch'
set -euo pipefail
EXP_ID="REANCHOR-CEILING-ANNEVO-M7"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ANNEVO_DIR="${ROOT}/refs/repos/annevo-2026"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; PRED_DIR="${OUT_DIR}/predictions"; METRIC_DIR="${OUT_DIR}/metrics"
LOG_DIR="${OUT_DIR}/logs"; STATUS_FILE="${OUT_DIR}/STATUS"
SUBSETS="${ROOT}/outputs/SCREENREF-tiberius_like-ho-s0/eval_subsets"   # deterministic TEST split (shared by all runs)
mkdir -p "${PRED_DIR}" "${METRIC_DIR}" "${LOG_DIR}"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate annevo; set -u
export TMPDIR="/tmp/annevo_${SLURM_JOB_ID:-$$}"; mkdir -p "${TMPDIR}"
cd "${ROOT}"
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())"

run_annevo() {
  local slug="$1" clade="$2" lineage="$3"
  local genome="data/m1_screen/${slug}/genome.fa"
  local pred="${PRED_DIR}/${slug}.gff"
  echo "[${EXP_ID}] ANNEVO clade=${clade} for ${slug}"
  ( cd "${ANNEVO_DIR}" && python annotation.py -g "${ROOT}/${genome}" \
      -m "saved_model/ANNEVO_${clade}.pt" -l "${lineage}" -o "${pred}" \
      --batch_size 32 -t "${SLURM_CPUS_PER_TASK:-16}" --show_log )
  # eval restricted to TEST chromosomes (eval_subsets genome) + new full-transcript ruler
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${SUBSETS}/${slug}/reference.gff3" \
    --prediction-gtf "${pred}" \
    --genome-fasta "${SUBSETS}/${slug}/genome.fa" \
    --output-json "${METRIC_DIR}/${slug}.metrics.json" \
    --experiment-id "${EXP_ID}_${slug}" --profile screen --span-mode cds
}

run_annevo "arabidopsis_thaliana" "Magnoliopsida" "Magnoliopsida"
run_annevo "gallus_gallus" "Aves" "Aves"

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" "${METRIC_DIR}/gallus_gallus.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
