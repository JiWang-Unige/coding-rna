#!/usr/bin/env bash
#SBATCH --job-name=M22GTVSMK
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --time=03:50:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu034,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
SEED=0
SP="arabidopsis_thaliana"
EXP_ID="M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-SMOKE-s${SEED}"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
STATUS_FILE="${OUT_DIR}/STATUS"

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
  echo "parent_exp_id=M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s0"
  echo "model_name=GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
  echo "objective=gb_tversky_smoke_noncrf"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
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
  --species "data/m1_screen/${SP}" \
  --exp-id "${EXP_ID}" --out-dir "${OUT_DIR}" \
  --model-name GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview --model-task token_classification \
  --window 6144 --sample-fraction 0.12 --epochs 1 --patience 1 \
  --batch-size 1 --eval-batch 4 --predict-batch 4 --lr 8e-4 --lora-lr 2e-5 --seed "${SEED}" \
  --class-weighting sqrt_inv --loss gb_tversky \
  --gb-tversky-lambda 2.0 --gb-tversky-alpha 0.8 --gb-tversky-beta 0.2 \
  --decoder none \
  --postproc constrained --min-cds-len 90 --max-fill-gap 20 \
  --lora-r 8 --lora-alpha 16 --lora-dropout 0.05 \
  --lora-targets q_proj,k_proj,v_proj,o_proj --attn-implementation sdpa --bf16 \
  --limit-train-windows 8 --limit-val-windows 4 --limit-test-seqids 1 --limit-predict-windows 8 \
  --pretokenize-windows --train-progress-every 1 --eval-progress-every 1 \
  --save-raw-scores

python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT_DIR}/eval_subsets/${SP}/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/${SP}.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/${SP}/genome.fa" \
  --output-json "${METRIC_DIR}/${SP}.metrics.json" \
  --experiment-id "${EXP_ID}_${SP}" --profile smoke --span-mode cds

python scripts/aggregate_gene_body_metrics.py \
  --metrics "${METRIC_DIR}/${SP}.metrics.json" \
  --output-json "${METRIC_DIR}/metrics.json" \
  --experiment-id "${EXP_ID}" --profile smoke

echo "COMPLETED" > "${STATUS_FILE}"

trap - EXIT
echo "DONE ${EXP_ID}"
