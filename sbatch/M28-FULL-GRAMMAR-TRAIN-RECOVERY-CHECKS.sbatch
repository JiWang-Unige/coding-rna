#!/usr/bin/env bash
#SBATCH --job-name=M28TRAINRECCHECK
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-FULL-GRAMMAR-TRAIN-RECOVERY-CHECKS/logs/%x_%j.out
#SBATCH --error=outputs/M28-FULL-GRAMMAR-TRAIN-RECOVERY-CHECKS/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OUT=outputs/M28-FULL-GRAMMAR-TRAIN-RECOVERY-CHECKS
test ! -e "$OUT/STATUS"
printf 'RUNNING\n' > "$OUT/STATUS"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$OUT/STATUS"; fi; exit "$rc"' EXIT
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/check_full_grammar_train_recovery.py
printf 'COMPLETED\n' > "$OUT/STATUS"
