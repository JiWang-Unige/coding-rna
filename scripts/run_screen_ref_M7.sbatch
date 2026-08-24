#!/usr/bin/env bash
#SBATCH --job-name=SREFM7
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=outputs/screen_ref_logs/%x_%j.out
#SBATCH --error=outputs/screen_ref_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# REANCHOR-HELDOUT-M7 CK3 — held-out same-budget anchor ref: random-init tiberius_like trained on
# {arabidopsis_thaliana (full), gallus_gallus (<=20Mb subset)} under the UNIFIED screen protocol,
# eval new full-transcript intergenic ruler (--span-mode cds), aggregate, validate. NON-CLAIM.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_screen_ref_M7.sbatch <seed>'
set -euo pipefail
SEED="${1:?usage: run_screen_ref_M7.sbatch <seed>}"
MODEL="tiberius_like"
EXP_ID="SCREENREF-${MODEL}-ho-s${SEED}"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"; SP_GG="data/m1_screen/gallus_gallus"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/screen_ref_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u
export TMPDIR="/tmp/srefm7_${SLURM_JOB_ID:-$$}"; mkdir -p "${TMPDIR}"
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())"

python -m src.screen_anchor.train_screen_ref \
  --model "${MODEL}" --decoder softmax --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --species "${SP_ARA}" "${SP_GG}" \
  --window 2048 --sample-fraction 0.3 --epochs 8 --patience 3 \
  --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv

for sp in arabidopsis_thaliana gallus_gallus; do
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" \
    --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" \
    --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" "${METRIC_DIR}/gallus_gallus.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
set +e
python scripts/validate_goal.py --goal ACTIVE_GOAL.json \
  --metrics "${METRIC_DIR}/metrics.json" --profile screen --run-status "${STATUS_FILE}" \
  > "${METRIC_DIR}/validate_goal.json"
echo "$?" > "${METRIC_DIR}/validate_goal.exit_code"
set -e
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
