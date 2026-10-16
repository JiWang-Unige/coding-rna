#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/M27-SINGLE-EDGE-R1/logs/%x_%j.out
#SBATCH --error=outputs/M27-SINGLE-EDGE-R1/logs/%x_%j.err
#SBATCH --job-name=M27EDGE
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00
set -euo pipefail
cd /home/users/j/jwang/coding-rna
M27_OUT="$(pwd -P)/outputs/M27-SINGLE-EDGE-R1"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate annevo
set -u
test "$(cat "$M27_OUT/PREP_STATUS")" = COMPLETED
test ! -e "$M27_OUT/STATUS"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$M27_OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$M27_OUT/STATUS"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=8 CUDA_VISIBLE_DEVICES=""
python -u scripts/experiments/M27-ALLELIC-INTEGRITY/edge_run.py --output-dir "$M27_OUT"
printf 'COMPLETED\n' > "$M27_OUT/STATUS"
