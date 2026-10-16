#!/usr/bin/env bash
#SBATCH --job-name=M28CODONENDPOINTCHECKS
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --no-requeue
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD:$PWD/tests:$PWD/scripts/experiments/M28-METHOD-RESTART"
out=outputs/M28-FULL-GRAMMAR-CODON-TRAIN-CHECKS/$SLURM_JOB_ID
mkdir -p "$out"
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$out/STATUS"; fi; exit "$rc"' EXIT
python -u -m unittest discover -s tests -p test_m28_full_grammar_codon_train.py -v > "$out/tests.log" 2>&1
printf 'COMPLETED\n' > "$out/STATUS"
