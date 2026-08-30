#!/usr/bin/env bash
#SBATCH --job-name=M25RSMOKE
#SBATCH --partition=private-teodoro-gpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --output=outputs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0-SMOKE/logs/%x_%j.out
#SBATCH --error=outputs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0-SMOKE/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/home/users/j/jwang/coding-rna"
EXP_ID="M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0-SMOKE"
CONFIG="${ROOT}/configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
STATUS_FILE="${OUT_DIR}/STATUS"

cd "${ROOT}"
mkdir -p "${OUT_DIR}/logs"
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
set +u
conda activate generanno
set -u
export HF_HOME="${HOME}/.cache/huggingface"
export PYTHONUNBUFFERED=1

python -u -m src.foundation_probe.train_generanno_structural_heads \
  --config "${CONFIG}" \
  --exp-id "${EXP_ID}" \
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa \
  --out-dir "${OUT_DIR}" \
  --setaria-fasta data/m25_setaria/GCF_005286985.2/genome.fa \
  --model-name GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview \
  --model-task token_classification \
  --window 6144 \
  --sample-fraction 0.12 \
  --epochs 3 \
  --batch-size 1 \
  --predict-batch 1 \
  --head-lr 8e-4 \
  --lora-lr 2e-5 \
  --seed 0 \
  --fp-lambda 2.5 \
  --lora-r 8 \
  --lora-alpha 16 \
  --lora-dropout 0.05 \
  --lora-targets q_proj,k_proj,v_proj,o_proj \
  --attn-implementation sdpa \
  --bf16 \
  --limit-train-windows 1536 \
  --limit-val-windows 768 \
  --development-smoke

echo "PASSED_SMOKE" > "${STATUS_FILE}"
trap - EXIT
