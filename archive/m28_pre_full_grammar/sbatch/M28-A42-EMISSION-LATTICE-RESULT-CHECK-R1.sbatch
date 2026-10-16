#!/usr/bin/env bash
#SBATCH --job-name=M28A42CHECKR1
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART" CUDA_VISIBLE_DEVICES=""
: "${POTENTIAL_JOB:?Actual original frozen potential job required}" "${DECODE_JOB:?Actual original complete 32-shard decode required}" "${READOUT_JOB:?Actual closed R1 TRAIN label readout required}"
out=outputs/M28-A42-EMISSION-LATTICE-RESULT-CHECK-R1
test ! -e "$out"
mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u -c 'import os; import torch; import run_emission_lattice_train as D; import verify_emission_lattice_train as V; torch.set_num_threads(2); D.OUT=D.ROOT/"outputs/M28-A42-EMISSION-LATTICE-TRAIN-READOUT-R1"; D.CHECK=D.ROOT/"outputs/M28-A42-EMISSION-LATTICE-RESULT-CHECK-R1"; V.main(os.environ["POTENTIAL_JOB"], os.environ["DECODE_JOB"], os.environ["READOUT_JOB"])'
printf 'COMPLETED\n' > "$out/STATUS"
