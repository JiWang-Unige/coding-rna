#!/usr/bin/env bash
#SBATCH --job-name=M28AQJOINT
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:01:36
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
# One allowed assertion repair: 8 joint tests; the 17 raw tests already passed.
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=""
out=outputs/M28-A-QUERY-JOINT-R2
test ! -e "$out"
mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"
printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python -B -m unittest discover -s tests -p test_m28_query_joint.py -v
PYTHONPATH=/home/users/j/jwang/coding-rna python -B -c 'import json; from src.m28.query_model import GeneQueryModel; m=GeneQueryModel(2048,16,128); print(json.dumps({"scope":"untrained parameter count only","feature_dim":2048,"num_queries":16,"hidden":128,"trainable_parameters":sum(p.numel() for p in m.parameters()),"optimizer_steps":0}))'
printf 'COMPLETED\n' > "$out/STATUS"
