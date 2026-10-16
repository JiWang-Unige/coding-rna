#!/usr/bin/env bash
#SBATCH --job-name=M28COMPLETECOREDRAWS
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
  CHECKS) OUT=outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-PREPARATION-CHECKS; args=(--fixture-only) ;;
  PREPARE) OUT=outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-PREPARATION; args=() ;;
  *) exit 2 ;;
esac
test ! -e "$OUT/STATUS"
printf 'RUNNING\n' > "$OUT/STATUS"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$OUT/STATUS"; fi; exit "$rc"' EXIT
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/prepare_full_grammar_complete_core_draws.py "${args[@]}"
printf 'COMPLETED\n' > "$OUT/STATUS"
