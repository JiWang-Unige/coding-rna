#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONALIBIO
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --nodelist=gpu035
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:02:00
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export CUDA_VISIBLE_DEVICES=""
unset PYTHONPATH
export OUT="$PWD/outputs/M28-CARBON-A-LIBRARY-IO/$SLURM_JOB_ID"
test ! -e "$OUT"; mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
/home/users/j/jwang/.conda/envs/generanno/bin/python -u scripts/experiments/M28-METHOD-RESTART/probe_carbon_a_library_io.py
printf 'COMPLETED\n' > "$OUT/STATUS"
