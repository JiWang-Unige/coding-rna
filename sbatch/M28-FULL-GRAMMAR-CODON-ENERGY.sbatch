#!/usr/bin/env bash
#SBATCH --job-name=M28CODONENERGY
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 TOKENIZERS_PARALLELISM=false
export HF_HOME="/home/users/j/jwang/.cache/huggingface" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
export PYTHONPATH="$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART:$PWD/tests"
case "$1" in
  checks) out=outputs/M28-FULL-GRAMMAR-CODON-ENERGY-CHECKS ;;
  task)
    seqids=(NC_003070.9 NC_003071.7 NC_003076.8 NC_089043.1 NC_089044.1 NC_089046.1)
    out="outputs/M28-FULL-GRAMMAR-CODON-ENERGY-DIAGNOSTIC/${seqids[$SLURM_ARRAY_TASK_ID]}" ;;
  aggregate) out=outputs/M28-FULL-GRAMMAR-CODON-ENERGY-DIAGNOSTIC ;;
  *) exit 2 ;;
esac
if [[ "$1" == aggregate ]]; then
  test ! -e "$out/summary.json" && test ! -e "$out/STATUS"
else
  test ! -e "$out"; mkdir -p "$out"
fi
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
if [[ "$1" == task ]]; then
  python -u scripts/experiments/M28-METHOD-RESTART/diagnose_full_grammar_codon_energy.py task --task "$SLURM_ARRAY_TASK_ID"
else
  python -u scripts/experiments/M28-METHOD-RESTART/diagnose_full_grammar_codon_energy.py "$1"
fi
printf 'COMPLETED\n' > "$out/STATUS"
