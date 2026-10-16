#!/usr/bin/env bash
#SBATCH --job-name=M28A30SAVEDLOSS
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=""
out=outputs/M28-A30-SAVED-SEARCH-LOSS-DIAGNOSTIC
test ! -e "$out"
trap 'rc=$?; if [[ $rc -ne 0 && -d "$out" ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
python -u scripts/experiments/M28-METHOD-RESTART/diagnose_saved_search_retention_records.py
