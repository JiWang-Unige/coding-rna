#!/usr/bin/env bash
#SBATCH --job-name=M25RGRAPHTEST
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.out
#SBATCH --error=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-RUN-GRAPH-R1
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED_TEST exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
trap 'exit 143' TERM INT
printf 'RUNNING_SYNTHETIC_TEST\n' > "$OUT/STATUS"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2
python -m pytest -q -p no:cacheprovider tests/test_m25r_run_graph.py --junitxml="$OUT/synthetic_tests.xml"
printf 'SYNTHETIC_TEST_PASS\n' > "$OUT/STATUS"
