#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONAPREQUAL
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
export SOURCE_NATIVE_JOB="${1:?pass the consumed native preflight job ID}"
export ARCHIVE_JOB="${2:?pass its completed software archive job ID}"
export PATH="$runtime/bin:$PATH" CUDA_VISIBLE_DEVICES=""
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
export OUT="$PWD/outputs/M28-CARBON-A-PREFLIGHT/$SLURM_JOB_ID"
test ! -e "$OUT"; mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
test "$( < "$PWD/outputs/M28-CARBON-A-SOFTWARE-ARCHIVE/$ARCHIVE_JOB/STATUS")" = COMPLETED
export TASK_TF_STAGE="$(mktemp -d /tmp/m28-carbon-a-tf.$SLURM_JOB_ID.XXXXXX)"
tar -xzf "$PWD/outputs/M28-CARBON-A-SOFTWARE-ARCHIVE/$ARCHIVE_JOB/transformers-4.56.2.tar.gz" -C "$TASK_TF_STAGE"
export PYTHONPATH="$TASK_TF_STAGE"
python -u scripts/experiments/M28-METHOD-RESTART/qualify_carbon_a_preflight.py
printf 'COMPLETED\n' > "$OUT/STATUS"
