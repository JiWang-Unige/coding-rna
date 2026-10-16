#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONARTQUAL
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
unset PYTHONPATH
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
export PATH="$runtime/bin:$PATH" OUT="$PWD/outputs/M28-CARBON-A-RUNTIME/$SLURM_JOB_ID"
test ! -e "$OUT"; mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -m pip check
python -m pip freeze > "$OUT/pip-freeze.txt"
cmp outputs/M28-CARBON-A-RUNTIME/13500426/pip-freeze.txt "$OUT/pip-freeze.txt"
python -u scripts/experiments/M28-METHOD-RESTART/qualify_carbon_a_runtime.py
printf 'COMPLETED\n' > "$OUT/STATUS"
