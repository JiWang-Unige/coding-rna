#!/usr/bin/env bash
#SBATCH --job-name=M28CODONCHECKS
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --no-requeue
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD:$PWD/tests"
OUT=outputs/M28-FULL-GRAMMAR-CODON-CHECKS/$SLURM_JOB_ID
mkdir -p "$OUT"
printf 'RUNNING\n' > "$OUT/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$OUT/STATUS"; fi; exit "$rc"' EXIT
python -u -m unittest discover -s tests -p test_m28_full_grammar_codon.py -v > "$OUT/tests.log" 2>&1
printf 'COMPLETED\n' > "$OUT/STATUS"
