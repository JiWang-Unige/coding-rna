#!/usr/bin/env bash
#SBATCH --job-name=SREFM8
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=outputs/screen_ref_logs/%x_%j.out
#SBATCH --error=outputs/screen_ref_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TB-GBF1-MULTICLASS-M8 CK3 — same-budget raw-DNA anchor (tiberius_like, 3-class) on the CLEAN
# held-out plants {arabidopsis, rice} = the gbF1 BAR (M8 candidate must not be < this on gbF1).
# Usage: ssh baobab 'cd .. && sbatch scripts/run_screen_ref_M8.sbatch <seed>'
set -euo pipefail
SEED="${1:?usage: <seed>}"; MODEL="tiberius_like"
EXP_ID="SCREENREF-${MODEL}-m8clean-s${SEED}"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"; SP_RICE="data/m1_screen/oryza_sativa"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/screen_ref_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u
export TMPDIR="/tmp/srefm8_${SLURM_JOB_ID:-$$}"; mkdir -p "${TMPDIR}"

python -m src.screen_anchor.train_screen_ref \
  --model "${MODEL}" --decoder softmax --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --species "${SP_ARA}" "${SP_RICE}" \
  --window 2048 --sample-fraction 0.3 --epochs 8 --patience 3 \
  --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv

for sp in arabidopsis_thaliana oryza_sativa; do
  python scripts/eval_gene_body_mask.py --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" "${METRIC_DIR}/oryza_sativa.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
