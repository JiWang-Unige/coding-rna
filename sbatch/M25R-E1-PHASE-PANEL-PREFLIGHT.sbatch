#!/usr/bin/env bash
#SBATCH --job-name=M25RPHASEPREF
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=outputs/M25R-E1-PHASE-GATE-PANEL/logs/%x_%j.out
#SBATCH --error=outputs/M25R-E1-PHASE-GATE-PANEL/logs/%x_%j.err

set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M25R-E1-PHASE-GATE-PANEL
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export CUDA_VISIBLE_DEVICES=""
export PYTHONDONTWRITEBYTECODE=1
python -m pytest -q -p no:cacheprovider tests/test_m25r_phase_panel_preflight.py
python -u scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/phase_panel_preflight.py
