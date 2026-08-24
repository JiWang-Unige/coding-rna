#!/usr/bin/env bash
#SBATCH --job-name=TACF5
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=03:55:00
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TA-COHERENCE-FIX-M5 — FP-FRAGFIX-CONSTR: FPLOSS (FP-aware loss) + constrained-decode post-proc
# (de-fragment) on frozen SegmentNT features, same-budget screen, NEW ruler. NON-CLAIM.
# Usage: sbatch [--partition .. --gres ..] scripts/run_TA-COHERENCE-FIX-M5.sbatch <seed> [epochs]
set -euo pipefail
SEED="${1:?usage: run_TA-COHERENCE-FIX-M5.sbatch <seed> [epochs]}"
EPOCHS="${2:-8}"
EXP_ID="FP-FRAGFIX-CONSTR-s${SEED}"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
CACHE="${ROOT}/outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_SAC="data/m1_screen/saccharomyces_cerevisiae"; SP_DMEL="data/m1_screen/drosophila_melanogaster"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u

python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP_SAC}" "${SP_DMEL}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --head convlstm --window 2048 --sample-fraction 0.3 --epochs "${EPOCHS}" --patience 3 \
  --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained

for sp in saccharomyces_cerevisiae drosophila_melanogaster; do
  python scripts/eval_gene_body_mask.py --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/saccharomyces_cerevisiae.metrics.json" "${METRIC_DIR}/drosophila_melanogaster.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
set +e
python scripts/validate_goal.py --goal ACTIVE_GOAL.json --metrics "${METRIC_DIR}/metrics.json" \
  --profile screen --run-status "${STATUS_FILE}" > "${METRIC_DIR}/validate_goal.json"
echo "$?" > "${METRIC_DIR}/validate_goal.exit_code"
set -e
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
