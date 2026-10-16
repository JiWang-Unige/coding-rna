#!/usr/bin/env bash
#SBATCH --job-name=M28CARBRAWREVALID
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
runtime="$PWD/refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
export PATH="$runtime/bin:$PATH" CUDA_VISIBLE_DEVICES=""
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
unset PYTHONPATH
export OUT="$PWD/outputs/M28-CARBON-A-TRAIN-RAW-VALIDATION/NC_003071.7"
test ! -e "$OUT"; mkdir -p "$OUT"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/revalidate_carbon_a_task1.py
printf 'COMPLETED\n' > "$OUT/STATUS"
