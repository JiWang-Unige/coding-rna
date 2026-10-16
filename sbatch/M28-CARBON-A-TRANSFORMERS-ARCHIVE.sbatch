#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONATFARCHIVE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
export PATH="$runtime/bin:$PATH"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
export OUT="$PWD/outputs/M28-CARBON-A-SOFTWARE-ARCHIVE/$SLURM_JOB_ID"
test ! -e "$OUT"
mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/prepare_carbon_a_transformers_archive.py
printf 'COMPLETED\n' > "$OUT/STATUS"
