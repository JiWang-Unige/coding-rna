#!/usr/bin/env bash
#SBATCH --job-name=M28A26FITREAD
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD"
fit_out=outputs/M28-A26-PHASE-EVIDENCE-FIT
readout_out=outputs/M28-A26-PHASE-EVIDENCE-TRAIN-READOUT
test ! -e "$fit_out"; test ! -e "$readout_out"
mkdir "$fit_out" "$readout_out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$readout_out/STATUS"; if [[ ! -s "$fit_out/summary.json" ]]; then printf "FAILED exit=%s\n" "$rc" > "$fit_out/STATUS"; fi; fi' EXIT
printf 'RUNNING\n' > "$fit_out/STATUS"
printf 'RUNNING\n' > "$readout_out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$fit_out/JOBID"
printf '%s\n' "$SLURM_JOB_ID" > "$readout_out/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/run_phase_evidence_chain_set.py --mode fit-readout
printf 'COMPLETED\n' > "$readout_out/STATUS"

