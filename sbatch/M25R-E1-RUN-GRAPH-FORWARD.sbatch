#!/usr/bin/env bash
#SBATCH --job-name=M25RGRAPHFWD
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.out
#SBATCH --error=outputs/M25R-E1-RUN-GRAPH-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-RUN-GRAPH-R1
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED_FORWARD_OR_A exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
trap 'exit 143' TERM INT
test "$(tr -d '\n' < "$OUT/STATUS")" = SYNTHETIC_TEST_PASS
test ! -e "$OUT/resolved_run.json"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export HF_HOME=/home/users/j/jwang/.cache/huggingface
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8
python -m pytest -q -p no:cacheprovider tests/test_m25r_run_graph_export.py --junitxml="$OUT/export_tests.xml"
printf 'RUNNING_INFERENCE\n' > "$OUT/STATUS"
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py infer
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py replay
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py evaluate --arm A
