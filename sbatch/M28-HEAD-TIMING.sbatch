#!/usr/bin/env bash
#SBATCH --job-name=M28HEADTIME
#SBATCH --no-requeue
#SBATCH --partition=shared-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu034,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-HEAD-TIMING-R3/logs/%x_%j.out
#SBATCH --error=outputs/M28-HEAD-TIMING-R3/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
out=outputs/M28-HEAD-TIMING-R3/result
mkdir -p "$out"
test ! -e "$out/summary.json"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -u scripts/experiments/M28-METHOD-RESTART/head_timing.py --feature-dir outputs/M28-GLM-FEATURE-SMOKE-R3/repair1 --output-dir "$out"
printf 'COMPLETED\n' > "$out/STATUS"
