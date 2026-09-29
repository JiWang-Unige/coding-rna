#!/usr/bin/env bash
#SBATCH --job-name=M28INTEGRATE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=outputs/M28-B1-INTEGRATION-R3/logs/%x_%j.out
#SBATCH --error=outputs/M28-B1-INTEGRATION-R3/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
species="$1"
out="outputs/M28-B1-INTEGRATION-R3/$species"
mkdir -p "$out"
test ! -e "$out/summary.json"
export NUMBA_CACHE_DIR="$PWD/$out/numba-cache"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -m unittest discover -s tests -p test_m28_candidates.py -v
python -u scripts/experiments/M28-METHOD-RESTART/integration_smoke.py --species "$species" --feature-dir outputs/M28-GLM-FEATURE-SMOKE-R3/repair1 --output-dir "$out"
printf 'COMPLETED\n' > "$out/STATUS"
