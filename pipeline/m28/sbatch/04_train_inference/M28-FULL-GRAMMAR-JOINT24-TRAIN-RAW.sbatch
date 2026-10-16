#!/usr/bin/env bash
#SBATCH --job-name=M28JOINT24TRAINRAW
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --array=0-5%3
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%A_%a.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%A_%a.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 TOKENIZERS_PARALLELISM=false
export HF_HOME="/home/users/j/jwang/.cache/huggingface" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
export PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts"
seqids=(NC_003070.9 NC_003071.7 NC_003076.8 NC_089043.1 NC_089044.1 NC_089046.1)
out="outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-RAW/${seqids[$SLURM_ARRAY_TASK_ID]}"
test ! -e "$out"; mkdir -p "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -u pipeline/m28/scripts/infer_full_grammar_joint24_train.py raw --task "$SLURM_ARRAY_TASK_ID" --fit-job "$1"
printf 'COMPLETED\n' > "$out/STATUS"
