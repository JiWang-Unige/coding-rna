#!/usr/bin/env bash
#SBATCH --job-name=M28FULLPRODRAW
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:02:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES="" PYTHONPATH="$PWD/pipeline/m28:$PWD/tests:$PWD/pipeline/m28/scripts"
out=outputs/M28-FULL-GRAMMAR-PRODUCTION-RAW-CHECKS
test ! -e "$out"; mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u pipeline/m28/scripts/check_full_grammar_production_raw.py
printf 'COMPLETED\n' > "$out/STATUS"
