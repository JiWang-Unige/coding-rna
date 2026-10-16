#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONARTVERIFY
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
unset PYTHONPATH
out=outputs/M28-CARBON-A-RUNTIME/$SLURM_JOB_ID
runtime=refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1
install=outputs/M28-CARBON-A-RUNTIME/13500092
test ! -e "$out"; mkdir -p "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
test -x "$runtime/bin/python"
test -s "$install/pip-freeze.txt"
export PATH="$PWD/$runtime/bin:$PATH" CARBON_A_INSTALL_JOB=13500092
# The installation and its pip check succeeded before the first job timed out.
# Reuse it unchanged; this one bounded repair performs verification only.
python -m pip check
python -m pip freeze > "$out/pip-freeze.txt"
cmp "$install/pip-freeze.txt" "$out/pip-freeze.txt"
python -u scripts/experiments/M28-METHOD-RESTART/verify_carbon_a_runtime.py "$runtime" "$out"
printf 'COMPLETED\n' > "$out/STATUS"
