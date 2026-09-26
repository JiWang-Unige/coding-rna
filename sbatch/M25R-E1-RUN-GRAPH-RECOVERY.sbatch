#!/usr/bin/env bash
#SBATCH --job-name=M25RGRAPHRECOVER
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY/logs/%x_%j.out
#SBATCH --error=outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY
SOURCE_RUN=outputs/M25R-E1-RUN-GRAPH-R1
SCRIPT=scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/run_graph_run.py
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED_CPU_RECOVERY exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
trap 'exit 143' TERM INT
test ! -e "$OUT/A_models.json"
test ! -e "$OUT/B_models.json"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8
export CUDA_VISIBLE_DEVICES=""
printf 'RUNNING_RECOVERY_TESTS\n' > "$OUT/STATUS"
python -m pytest -q -p no:cacheprovider tests/test_m25r_run_graph_recovery.py tests/test_m25r_run_graph_export.py --junitxml="$OUT/recovery_tests.xml"
test "$(jq -r .windows_forwarded "$SOURCE_RUN/inference_complete.json")" = 17384
cp -- "$SOURCE_RUN/A_models.json" "$SOURCE_RUN/A_predictions.gff3" "$SOURCE_RUN/A_prediction_complete.json" "$OUT/"
printf 'RUNNING_A_EVALUATION\n' > "$OUT/STATUS"
python -u "$SCRIPT" evaluate --arm A --output-dir "$OUT"
printf 'RUNNING_B_DECODE\n' > "$OUT/STATUS"
python -u "$SCRIPT" decode --output-dir "$OUT" --score-source-run "$SOURCE_RUN"
printf 'RUNNING_B_EVALUATION\n' > "$OUT/STATUS"
python -u "$SCRIPT" evaluate --arm B --output-dir "$OUT"
