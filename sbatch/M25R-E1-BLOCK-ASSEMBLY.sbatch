#!/usr/bin/env bash
#SBATCH --job-name=M25RBLKCPU
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/M25RBLKCPU_%j.out
#SBATCH --error=logs/M25RBLKCPU_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-BLOCK-ASSEMBLY-REACHABILITY
mkdir "$OUT"
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/block_assembly_reachability.py
