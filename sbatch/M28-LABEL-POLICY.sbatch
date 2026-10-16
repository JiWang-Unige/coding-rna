#!/usr/bin/env bash
#SBATCH --job-name=M28LABELR4
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=outputs/M28-LABEL-POLICY-R4/logs/%x_%j.out
#SBATCH --error=outputs/M28-LABEL-POLICY-R4/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
species="$1"
out="outputs/M28-LABEL-POLICY-R4/$species"
mkdir -p "$out"
test ! -e "$out/summary.json"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -m pytest -q tests/test_m28_common_ruler.py tests/test_m28_labels.py tests/test_m28_label_alignment.py
python -u scripts/experiments/M28-METHOD-RESTART/align_labels.py --species "$species" --output-dir "$out"
printf 'COMPLETED\n' > "$out/STATUS"
