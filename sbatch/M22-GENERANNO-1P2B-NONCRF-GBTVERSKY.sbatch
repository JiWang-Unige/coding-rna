#!/usr/bin/env bash
#SBATCH --job-name=M22GENGTV
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --time=47:50:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/fp_segnt_logs/%x_%A_%a.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%A_%a.err
#SBATCH --exclude=gpu023,gpu024,gpu034,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
SEED="${SLURM_ARRAY_TASK_ID:-0}"
EXP_ID="M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s${SEED}"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
STATUS_FILE="${OUT_DIR}/STATUS"
SP_LIST=(arabidopsis_thaliana oryza_sativa)

mkdir -p "${METRIC_DIR}" "${ROOT}/outputs/fp_segnt_logs"
python "${ROOT}/scripts/pre_submit_gate.py" --exp-id "${EXP_ID}" "${ROOT}" --format text

if [[ -e "${STATUS_FILE}" ]] && grep -qx "COMPLETED" "${STATUS_FILE}"; then
  echo "${EXP_ID} already COMPLETED; refusing to overwrite."
  exit 0
fi
if [[ -d "${OUT_DIR}" ]] && find "${OUT_DIR}" -mindepth 1 -type f ! -name code_review_gate.json ! -name STATUS -print -quit | grep -q .; then
  echo "Refusing to overwrite non-empty output dir: ${OUT_DIR}" >&2
  exit 3
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
set +u
conda activate generanno
set -u
export HF_HOME="${HOME}/.cache/huggingface"
export PYTHONUNBUFFERED=1

{
  echo "exp_id=${EXP_ID}"
  echo "parent_exp_id=M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s${SEED}"
  echo "model_name=GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
  echo "objective=gb_tversky_lambda2_alpha0.8_beta0.2_noncrf_raw_scores"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
  echo "host=$(hostname)"
  date -Is
  python -V
  python - <<'PY'
import torch, transformers, peft
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("transformers", transformers.__version__, "peft", peft.__version__)
PY
} > "${OUT_DIR}/env.txt" 2>&1

python -u -m src.foundation_probe.train_generanno_lora_3class \
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --model-name GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview --model-task token_classification \
  --window 6144 --sample-fraction 0.12 --epochs 3 --patience 1 \
  --batch-size 1 --eval-batch 16 --predict-batch 16 --lr 8e-4 --lora-lr 2e-5 --seed "${SEED}" \
  --class-weighting sqrt_inv --loss gb_tversky \
  --gb-tversky-lambda 2.0 --gb-tversky-alpha 0.8 --gb-tversky-beta 0.2 \
  --decoder none \
  --postproc constrained --min-cds-len 90 --max-fill-gap 20 \
  --lora-r 8 --lora-alpha 16 --lora-dropout 0.05 \
  --lora-targets q_proj,k_proj,v_proj,o_proj --attn-implementation sdpa --bf16 \
  --limit-train-windows 1536 --limit-val-windows 768 \
  --pretokenize-windows --train-progress-every 200 --eval-progress-every 20 \
  --save-raw-scores

metric_paths=()
for SP in "${SP_LIST[@]}"; do
  python scripts/eval_gene_body_mask.py \
    --reference-gtf "${OUT_DIR}/eval_subsets/${SP}/reference.gff3" \
    --prediction-gtf "${OUT_DIR}/predictions/${SP}.gff" \
    --genome-fasta "${OUT_DIR}/eval_subsets/${SP}/genome.fa" \
    --output-json "${METRIC_DIR}/${SP}.metrics.json" \
    --experiment-id "${EXP_ID}_${SP}" --profile screen --span-mode cds
  metric_paths+=("${METRIC_DIR}/${SP}.metrics.json")
done

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${metric_paths[@]}" \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" --profile screen

python - <<'PY' "${METRIC_DIR}/metrics.json" "${METRIC_DIR}/m22_promotion_gate.json"
import json
import sys

metrics_path, out_path = sys.argv[1:3]
with open(metrics_path) as handle:
    metrics = json.load(handle)

gbf1 = float(metrics.get("gene_body_F1_unconstrained", float("nan")))
fpr = float(metrics.get("intergenic_FPR", float("nan")))
gcount = float(metrics.get("predicted_gene_count_ratio_vs_reference", float("nan")))
checks = {
    "hard_fpr_le_0p01": fpr <= 0.01,
    "gbf1_gt_m19_s1": gbf1 > 0.8815481561599101,
    "gene_count_le_1p25": gcount <= 1.25,
}
payload = {
    "exp_id": metrics.get("experiment_id"),
    "purpose": "M22-specific continue/abandon gate; stricter than generic screen validate_goal",
    "reference": {
        "m19_s1_gene_body_F1_unconstrained": 0.8815481561599101,
        "hard_intergenic_FPR_max": 0.01,
        "gene_count_ratio_max": 1.25,
    },
    "observed": {
        "gene_body_F1_unconstrained": gbf1,
        "intergenic_FPR": fpr,
        "predicted_gene_count_ratio_vs_reference": gcount,
    },
    "checks": checks,
    "promote": all(checks.values()),
}
with open(out_path, "w") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
print(json.dumps(payload, sort_keys=True))
PY

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
