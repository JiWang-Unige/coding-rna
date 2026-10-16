#!/usr/bin/env bash
#SBATCH --job-name=M28A31LOCALEVAL
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
generator=${1:?A3 or A23 required}
case "$generator" in A3|A23) ;; *) exit 2 ;; esac
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD" CUDA_VISIBLE_DEVICES=""
infer="outputs/M28-A31-${generator}-LOCAL-PATH-DEV-INFER"
out="outputs/M28-A31-${generator}-LOCAL-PATH-DEV-EVAL"
test ! -e "$out"
trap 'rc=$?; if [[ $rc -ne 0 && -d "$out" ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; printf "%s\n" "$SLURM_JOB_ID" > "$out/JOBID"; fi' EXIT
python -u scripts/experiments/M28-METHOD-RESTART/evaluate_local_path_caller.py \
  --generator "$generator" --inference-dir "$PWD/$infer" --output-dir "$PWD/$out"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
printf 'COMPLETED\n' > "$out/STATUS"
