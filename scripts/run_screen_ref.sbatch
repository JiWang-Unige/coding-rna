#!/usr/bin/env bash
#SBATCH --job-name=SCREENREF
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/screen_ref_logs/%x_%j.out
#SBATCH --error=outputs/screen_ref_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# M1-SAMEBUDGET-SCREEN-ANCHOR — train ONE random-init reference (model,seed) under the unified
# small-sample screen protocol, then eval per species with --span-mode cds + aggregate + validate.
# Usage: sbatch [--partition=.. --gres=..] scripts/run_screen_ref.sbatch <model> <seed>
#   model in {tiberius_like, helixer_like}
set -euo pipefail

MODEL="${1:?usage: run_screen_ref.sbatch <model> <seed> [decoder]}"
SEED="${2:?usage: run_screen_ref.sbatch <model> <seed> [decoder]}"
DECODER="${3:-softmax}"
if [[ "${DECODER}" == "softmax" ]]; then
  EXP_ID="SCREENREF-${MODEL}-s${SEED}"
else
  EXP_ID="SCREENREF-${MODEL}-${DECODER}-s${SEED}"
fi
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
STATUS_FILE="${OUT_DIR}/STATUS"
SP_SAC="data/m1_screen/saccharomyces_cerevisiae"
SP_DMEL="data/m1_screen/drosophila_melanogaster"

mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/screen_ref_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT

cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u                      # coding-rna now has mkl (via torch); its activate.d references unbound MKL_* vars
conda activate coding-rna
set -u
export TMPDIR="/tmp/screenref_${SLURM_JOB_ID:-$$}"; mkdir -p "${TMPDIR}"

python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())"

# ---- unified small-sample screen protocol (shared by refs AND future Track A candidates) ----
python -m src.screen_anchor.train_screen_ref \
  --model "${MODEL}" --decoder "${DECODER}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --species "${SP_SAC}" "${SP_DMEL}" \
  --window 2048 --sample-fraction 0.3 --epochs 8 --patience 3 \
  --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv

# ---- eval each species under the frozen CDS-span evaluator, then aggregate ----
for sp in saccharomyces_cerevisiae drosophila_melanogaster; do
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" \
    --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" \
    --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/saccharomyces_cerevisiae.metrics.json" \
            "${METRIC_DIR}/drosophila_melanogaster.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" --profile screen

set +e
python scripts/validate_goal.py --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" --profile screen --run-status "${STATUS_FILE}" \
  > "${METRIC_DIR}/validate_goal.json"
echo "$?" > "${METRIC_DIR}/validate_goal.exit_code"
set -e

echo "COMPLETED" > "${STATUS_FILE}"
trap - EXIT
