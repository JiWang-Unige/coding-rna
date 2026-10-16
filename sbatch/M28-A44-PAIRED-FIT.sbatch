#!/usr/bin/env bash
#SBATCH --job-name=M28A44FIT
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART"
arm=${1:?fixed A44 arm required}
checks_job=${2:?actual position-unit checks job required}
case "$arm" in
  A44_LOCAL) out=outputs/M28-A44-LOCAL-FIT ;;
  A44_JOINT) out=outputs/M28-A44-JOINT-FIT ;;
  *) printf 'Only A44_LOCAL/A44_JOINT allowed\n' >&2; exit 2 ;;
esac
test ! -e "$out"; mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u scripts/experiments/M28-METHOD-RESTART/fit_local_structured_emission.py --arm "$arm" --checks-job "$checks_job"
printf 'COMPLETED\n' > "$out/STATUS"
