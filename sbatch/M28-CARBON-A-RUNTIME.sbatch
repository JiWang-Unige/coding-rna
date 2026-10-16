#!/usr/bin/env bash
#SBATCH --job-name=M28CARBONART
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
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
assets=refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620
test ! -e "$out"; mkdir -p "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
test ! -e "$runtime"
test -f "$assets/requirements.txt"
export PIP_CACHE_DIR="$PWD/refs/carbon-a/runtimes/pip-cache"
export TMPDIR="$PWD/$out/tmp"
mkdir -p "$TMPDIR"
# Venv isolation: only the existing Python interpreter is reused, not its packages.
/home/users/j/jwang/.conda/envs/generanno/bin/python -m venv "$runtime"
export PATH="$PWD/$runtime/bin:$PATH"
python -m pip install --upgrade pip wheel setuptools
python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126
python -m pip install --only-binary=:all: transformers==4.56.2 -r "$assets/requirements.txt"
abi=$(python -c 'import torch; print("TRUE" if torch._C._GLIBCXX_USE_CXX11_ABI else "FALSE")')
wheel="flash_attn-2.8.3%2Bcu12torch2.7cxx11abi${abi}-cp310-cp310-linux_x86_64.whl"
python -m pip install "https://github.com/Dao-AILab/flash-attention/releases/download/v2.8.3/$wheel"
python -m pip check
python -m pip freeze > "$out/pip-freeze.txt"
python -u scripts/experiments/M28-METHOD-RESTART/verify_carbon_a_runtime.py "$runtime" "$out"
printf 'COMPLETED\n' > "$out/STATUS"
