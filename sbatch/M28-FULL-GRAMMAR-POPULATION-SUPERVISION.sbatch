#!/usr/bin/env bash
#SBATCH --job-name=M28POPULATIONFINE
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --no-requeue
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$PWD:$PWD/scripts:$PWD/scripts/experiments/M28-METHOD-RESTART"
case "${1:-}" in
  CHECKS) OUT=outputs/M28-FULL-GRAMMAR-POPULATION-SUPERVISION-CHECKS; args=(--fixture-only) ;;
  CENSUS) OUT=outputs/M28-FULL-GRAMMAR-POPULATION-SUPERVISION; args=() ;;
  *) exit 2 ;;
esac
test ! -e "$OUT/STATUS"
printf 'RUNNING\n' > "$OUT/STATUS"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$OUT/STATUS"; fi; exit "$rc"' EXIT
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/census_full_grammar_population_supervision.py "${args[@]}"
printf 'COMPLETED\n' > "$OUT/STATUS"
