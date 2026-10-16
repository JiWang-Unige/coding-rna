#!/usr/bin/env bash
#SBATCH --job-name=M28A36DEV
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:30:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD"
out=outputs/M28-A36-STRUCTURED-EMISSION-DEV-INFER
test ! -e "$out"
# Prepared only: actual TRAIN five-condition PASS and a separate DEV contract
# must precede submission. The Python entry verifies those source artifacts
# before model/FASTA access; never pre-create its absent output directory.
trap 'rc=$?; if [[ $rc -ne 0 && -d "$out" ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; printf "%s\n" "$SLURM_JOB_ID" > "$out/JOBID"; fi' EXIT
python -u scripts/experiments/M28-METHOD-RESTART/infer_structured_emission.py \
  --dev --output-dir "$PWD/$out"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
printf 'COMPLETED\n' > "$out/STATUS"
