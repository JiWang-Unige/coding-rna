#!/usr/bin/env bash
#SBATCH --job-name=M28A42DECODEPRE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --array=0-7%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:03:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%A_%a.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%A_%a.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART" CUDA_VISIBLE_DEVICES=""
: "${POTENTIAL_JOB:?actual potential job required}"
out=outputs/M28-A42-EMISSION-LATTICE-TRAIN-DECODE-PREFLIGHT/case_${SLURM_ARRAY_TASK_ID}
test ! -e "$out"; mkdir -p "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s_%s\n' "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" > "$out/JOBID"
python scripts/experiments/M28-METHOD-RESTART/preflight_emission_lattice_train.py --mode decode --case "$SLURM_ARRAY_TASK_ID" --potential-job "$POTENTIAL_JOB"
printf 'COMPLETED\n' > "$out/STATUS"
