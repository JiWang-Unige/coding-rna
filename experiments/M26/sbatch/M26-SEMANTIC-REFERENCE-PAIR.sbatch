#!/usr/bin/env bash
#SBATCH --job-name=M26SEMANTICPAIR
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=outputs/M26-REFERENCE-PAIR-R7-REPAIR1/logs/%x_%j.out
#SBATCH --error=outputs/M26-REFERENCE-PAIR-R7-REPAIR1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
OUT=outputs/M26-REFERENCE-PAIR-R7-REPAIR1
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
python -m pytest -q -p no:cacheprovider tests/test_m26_semantic_reference_pair.py --junitxml="$OUT/tests.xml"
python -u scripts/experiments/M26-SAME-SCOPE-MECHANISM/semantic_reference_pair.py --output-dir "$OUT"
printf 'COMPLETED\n' > "$OUT/STATUS"
