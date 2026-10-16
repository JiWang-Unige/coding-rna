#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONACTX24RAW
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --array=0-5%6
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%A_%a.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%A_%a.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
preflight_job=13503202
export ARCHIVE_JOB=13501094
source_archive="$PWD/outputs/M28-CARBON-A-SOFTWARE-ARCHIVE/$ARCHIVE_JOB"
export PATH="$runtime/bin:/opt/ebsofts/CUDA/12.6.0/bin:/opt/ebsofts/GCCcore/12.3.0/bin:$PATH"
export CUDA_HOME=/opt/ebsofts/CUDA/12.6.0 CXX=/opt/ebsofts/GCCcore/12.3.0/bin/g++
export TORCH_EXTENSIONS_DIR="$runtime/extensions" MAX_JOBS=2 TORCH_CUDA_ARCH_LIST=8.6
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export HF_HOME=/home/users/j/jwang/.cache/huggingface
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 TOKENIZERS_PARALLELISM=false
test "$(< "$PWD/outputs/M28-CARBON-A-PREFLIGHT/$preflight_job/STATUS")" = COMPLETED
test "$(< "$source_archive/STATUS")" = COMPLETED
export TASK_TF_STAGE="$(mktemp -d /tmp/m28-carbon-ctx24-tf.$SLURM_JOB_ID.XXXXXX)"
tar -xzf "$source_archive/transformers-4.56.2.tar.gz" -C "$TASK_TF_STAGE"
export PYTHONPATH="$TASK_TF_STAGE"
python -u scripts/experiments/M28-METHOD-RESTART/run_carbon_a_train.py \
    --task "$SLURM_ARRAY_TASK_ID" --runtime-job 13500755 --preflight-job "$preflight_job" --context-bp 24576
