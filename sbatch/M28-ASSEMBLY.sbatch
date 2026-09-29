#!/usr/bin/env bash
#SBATCH --job-name=M28ASSEMBLY
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=outputs/M28-ASSEMBLY-R4/logs/%x_%j.out
#SBATCH --error=outputs/M28-ASSEMBLY-R4/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
out="${1:-outputs/M28-ASSEMBLY-R4/result}"
mkdir -p "$out"
test ! -e "$out/STATUS"
export NUMBA_CACHE_DIR="$PWD/$out/numba-cache"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -m pytest -q tests/test_m28_assembly.py tests/test_m28_core.py tests/test_m28_c0_decode.py tests/test_m28_c0_gain.py
printf 'COMPLETED\n' > "$out/STATUS"
