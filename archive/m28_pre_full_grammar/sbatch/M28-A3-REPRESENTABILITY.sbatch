#!/usr/bin/env bash
#SBATCH --job-name=M28A3REPR
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:03:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONPATH="$PWD" PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
out=outputs/M28-A3-REPRESENTABILITY
test ! -e "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then mkdir -p "$out"; printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; printf "%s\n" "$SLURM_JOB_ID" > "$out/JOBID"; fi' EXIT
python -B -m unittest discover -s tests -p test_m28_representability.py -v
python -u scripts/experiments/M28-METHOD-RESTART/diagnose_representability.py
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
printf 'COMPLETED\n' > "$out/STATUS"
