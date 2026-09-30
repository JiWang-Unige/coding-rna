#!/usr/bin/env bash
#SBATCH --job-name=M28R6TRACE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
out="outputs/M28-R6-CANDIDATE-TRACE"
test ! -e "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "outputs/M28-PILOT-R5/logs/R6_FAILED_$SLURM_JOB_ID"; fi' EXIT
python -u scripts/experiments/M28-METHOD-RESTART/r6_trace.py --sample reports/M28-METHOD-RESTART/r6_trace_sample.json --output-dir "$out"
printf 'COMPLETED\n' > "$out/STATUS"
