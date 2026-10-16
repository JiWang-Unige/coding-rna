#!/usr/bin/env bash
#SBATCH --job-name=M28JOINTCEIL
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:01:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --no-requeue
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export CUDA_VISIBLE_DEVICES="" PYTHONDONTWRITEBYTECODE=1
out=outputs/M28-FULL-GRAMMAR-JOINT-REFERENCE-CEILING
test ! -e "$out/STATUS"
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED\n" > "$out/STATUS"; fi; exit "$rc"' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/diagnose_full_grammar_joint_reference_ceiling.py
printf 'COMPLETED\n' > "$out/STATUS"
