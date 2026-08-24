#!/usr/bin/env bash
#SBATCH --job-name=M8CAND
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# TB-GBF1-MULTICLASS-M8 CK3 — candidate on the CLEAN held-out plants {arabidopsis, rice}.
# scheme=mc  -> multi-class (8-class strand-aware) + CRF 8x8 transitions (the M8 gbF1-recovery bet)
# scheme=3c  -> 3-class FP-FRAGFIX-CONSTR config (the M7 baseline, re-run on the SAME clean species)
# both: fp_aware loss + constrained post-proc mfg=20/mcl=60, --save-raw-pred. eval new ruler.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_M8_cand.sbatch <mc|3c> <seed> [epochs]'
set -euo pipefail
SCHEME="${1:?usage: <mc|3c> <seed> [epochs]}"; SEED="${2:?seed}"; EPOCHS="${3:-8}"
if [[ "${SCHEME}" == "mc" ]]; then
  EXP_ID="M8-MC-CAND-s${SEED}"; LS="multiclass"; DEC="crf"
else
  EXP_ID="M8-3C-CAND-s${SEED}"; LS="3class"; DEC="none"
fi
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
CACHE="${ROOT}/outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"; METRIC_DIR="${OUT_DIR}/metrics"; STATUS_FILE="${OUT_DIR}/STATUS"
SP_ARA="data/m1_screen/arabidopsis_thaliana"; SP_RICE="data/m1_screen/oryza_sativa"
mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS_FILE}"; exit $rc; }
trap mark_failed EXIT
cd "${ROOT}"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u

python -m src.foundation_probe.train_probe_head \
  --cache "${CACHE}" --species "${SP_ARA}" "${SP_RICE}" --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --head convlstm --label-scheme "${LS}" --decoder "${DEC}" --window 2048 --sample-fraction 0.3 \
  --epochs "${EPOCHS}" --patience 3 --batch-size 16 --lr 1e-3 --seed "${SEED}" --class-weighting sqrt_inv \
  --loss fp_aware --fp-lambda 1.0 --postproc constrained --min-cds-len 60 --max-fill-gap 20 --save-raw-pred 2>&1 \
  | grep -ivE 'rematerialization|bfc_allocator|Current allocation'

for sp in arabidopsis_thaliana oryza_sativa; do
  python scripts/eval_gene_body_mask.py --reference-gtf "${OUT_DIR}/eval_subsets/${sp}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${sp}.gff" --genome-fasta "${OUT_DIR}/eval_subsets/${sp}/genome.fa" \
    --output-json "${METRIC_DIR}/${sp}.metrics.json" --experiment-id "${EXP_ID}_${sp}" --profile screen --span-mode cds
done
python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/arabidopsis_thaliana.metrics.json" "${METRIC_DIR}/oryza_sativa.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" --experiment-id "${EXP_ID}" --profile screen
echo "COMPLETED" > "${STATUS_FILE}"; trap - EXIT; echo "DONE ${EXP_ID}"
