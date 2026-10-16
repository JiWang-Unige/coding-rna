#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONALOCALTF
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
export ARCHIVE_JOB="${1:?pass the completed software-archive job ID}"
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
source_archive="$PWD/outputs/M28-CARBON-A-SOFTWARE-ARCHIVE/$ARCHIVE_JOB"
export OUT="$PWD/outputs/M28-CARBON-A-LOCAL-TF-PROBE/$SLURM_JOB_ID"
test ! -e "$OUT"
mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
test "$(< "$source_archive/STATUS")" = COMPLETED
export TASK_TF_STAGE="$(mktemp -d /tmp/m28-carbon-a-tf.$SLURM_JOB_ID.XXXXXX)"
tar -xzf "$source_archive/transformers-4.56.2.tar.gz" -C "$TASK_TF_STAGE"
export PATH="$runtime/bin:$PATH" PYTHONPATH="$TASK_TF_STAGE"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1
python -u scripts/experiments/M28-METHOD-RESTART/probe_carbon_a_local_transformers.py
printf 'COMPLETED\n' > "$OUT/STATUS"
