#!/usr/bin/env bash
#SBATCH --job-name=M28GRUEQ
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --exclude=gpu034
#SBATCH --output=outputs/M28-PACKED-GRU-CPU-R1/logs/%x_%j.out
#SBATCH --error=outputs/M28-PACKED-GRU-CPU-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
out=outputs/M28-PACKED-GRU-CPU-R1/result
mkdir -p "$out"
test ! -e "$out/STATUS"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
python -m pytest -q -s tests/test_m28_packed_chain_gru.py --junitxml="$out/junit.xml"
printf 'COMPLETED\n' > "$out/STATUS"
