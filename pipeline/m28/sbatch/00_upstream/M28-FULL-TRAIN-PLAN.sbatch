#!/usr/bin/env bash
#SBATCH --job-name=M28FTPLANR2
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%A_%a.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%A_%a.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
export PYTHONPATH="$PWD/pipeline/m28"
export CUDA_VISIBLE_DEVICES=""
out=outputs/M28-FULL-TRAIN-PLAN-R2
test ! -e "$out"
mkdir -p "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID:+_$SLURM_ARRAY_TASK_ID}" > "$out/JOBID"
python -B -m unittest discover -s tests -p test_m28_full_train.py -v
python -u pipeline/m28/scripts/prepare_full_train.py
printf 'COMPLETED\n' > "$out/STATUS"
