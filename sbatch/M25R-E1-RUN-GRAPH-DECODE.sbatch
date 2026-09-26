#!/usr/bin/env bash
#SBATCH --job-name=M25RGRAPHDECODE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.out
#SBATCH --error=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-RUN-GRAPH-R1
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED_DECODE_OR_EVALUATION exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
trap 'exit 143' TERM INT
test "$(tr -d '\n' < "$OUT/STATUS")" = A_FULL_REPLAY_PASS
test ! -e "$OUT/B_models.json"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8
printf 'RUNNING_DECODE\n' > "$OUT/STATUS"
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py decode
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py evaluate --arm B
