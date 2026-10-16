#!/usr/bin/env bash
#SBATCH --job-name=M28COUNTS
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=outputs/M28-METHOD-RESTART-R1/logs/%x_%j.out
#SBATCH --error=outputs/M28-METHOD-RESTART-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
python -u scripts/experiments/M28-METHOD-RESTART/verify_legacy_counts.py
