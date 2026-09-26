#!/usr/bin/env bash
#SBATCH --job-name=M26ENDPOINT
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=outputs/M26-ENDPOINT-BUDGET-R8/logs/%x_%j.out
#SBATCH --error=outputs/M26-ENDPOINT-BUDGET-R8/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 CUDA_VISIBLE_DEVICES=""
OUT=outputs/M26-ENDPOINT-BUDGET-R8
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
python -m pytest -q -p no:cacheprovider tests/test_m26_endpoint_budget.py --junitxml="$OUT/tests.xml"
python -u scripts/experiments/M26-SAME-SCOPE-MECHANISM/endpoint_budget.py --output-dir "$OUT"
printf 'COMPLETED\n' > "$OUT/STATUS"
