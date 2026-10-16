#!/usr/bin/env bash
#SBATCH --job-name=M28NATEVCHK
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:03:00
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --no-requeue
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=outputs/M28-NATIVE-EVALUATION-CHECKS
test ! -e "$OUT"
mkdir "$OUT"
printf 'RUNNING\n' > "$OUT/STATUS"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$OUT/STATUS"; fi; exit "$rc"' EXIT
printf '%s\n' "$SLURM_JOB_ID" > "$OUT/JOBID"
python -m unittest discover -s tests -p test_m28_native_baseline_evaluation.py -v
printf 'COMPLETED\n' > "$OUT/STATUS"
