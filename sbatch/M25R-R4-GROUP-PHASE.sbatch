#!/usr/bin/env bash
#SBATCH --job-name=M25RGROUPPHASE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE/logs/%x_%j.out
#SBATCH --error=outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE/logs/%x_%j.err

set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
python -m pytest -q -p no:cacheprovider tests/test_m25r_group_phase.py tests/test_m25r_case_correction.py tests/test_m25r_redecode_error_decomposition.py
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/group_phase_retrospective.py
