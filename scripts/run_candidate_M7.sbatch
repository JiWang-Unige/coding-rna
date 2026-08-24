#!/usr/bin/env bash
#SBATCH --job-name=CANDM7
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# REANCHOR-HELDOUT-M7 CK5 — held-out re-test of promote-ready candidate FP-FRAGFIX-CONSTR
# (frozen SegmentNT feats + FP-aware loss convlstm head + constrained post-proc mfg=20/mcl=60),
# IDENTICAL config, only species changed to {arabidopsis_thaliana, gallus_gallus(<=20Mb subset)}.
# --save-raw-pred so the constrained band can be VAL-re-selected on held-out if gene_count is off
# (VAL-only, TEST-once; no leakage). eval new full-transcript ruler. NON-CLAIM.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_candidate_M7.sbatch <seed> [epochs]'
set -euo pipefail
SEED="${1:?usage: <seed> [epochs]}"
EPOCHS="${2:-8}"
EXP_ID="FP-FRAGFIX-CONSTR-ho-s${SEED}"
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
CACHE="${ROOT}/outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"; SP_GG="data/m1_screen/gallus_gallus"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u

python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP_ARA}" "${SP_GG}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --head convlstm --window 2048 --sample-fraction 0.3 --epochs "${EPOCHS}" --patience 3 \
  --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 --save-raw-pred

for sp in arabidopsis_thaliana gallus_gallus; do
  python scripts/eval_gene_body_mask.py --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" "${METRIC_DIR}/gallus_gallus.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
