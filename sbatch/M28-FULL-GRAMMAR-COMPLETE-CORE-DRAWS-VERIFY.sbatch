#!/usr/bin/env bash
#SBATCH --job-name=M28CCVERIFY
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=""
export PYTHONPATH="$PWD:$PWD/scripts:$PWD/scripts/experiments/M28-METHOD-RESTART"
out=outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-VERIFICATION
test ! -e "$out/STATUS"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/fit_full_grammar_complete_core_draws.py verify --source-job "$1"
printf 'COMPLETED\n' > "$out/STATUS"
