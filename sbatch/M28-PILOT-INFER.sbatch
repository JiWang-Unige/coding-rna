#!/usr/bin/env bash
#SBATCH --job-name=M28INFER
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --exclude=gpu034
#SBATCH --output=outputs/M28-PILOT-R4/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R4/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
arm="$1"
out="outputs/M28-PILOT-R4/infer_$arm"
mkdir -p "$out"
test ! -e "$out/summary.json"
export NUMBA_CACHE_DIR="$PWD/$out/numba-cache"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -u scripts/experiments/M28-METHOD-RESTART/infer_pilot.py --arm "$arm" --pilot-dir outputs/M28-PILOT-R4 --output-dir "$out"
printf 'COMPLETED\n' > "$out/STATUS"
