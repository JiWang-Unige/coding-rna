#!/usr/bin/env bash
#SBATCH --job-name=M27PROBE
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=outputs/M27-ALLELIC-PREP-R2/logs/%x_%j.out
#SBATCH --error=outputs/M27-ALLELIC-PREP-R2/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
OUT=outputs/M27-ALLELIC-PREP-R2
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
apptainer exec --cleanenv --env PYTHONDONTWRITEBYTECODE=1 --bind "$(pwd -P):/work" --pwd /work \
    refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif \
    python scripts/experiments/M27-ALLELIC-INTEGRITY/runtime_probe.py "$OUT/tiberius_runtime.json"
printf 'COMPLETED\n' > "$OUT/STATUS"
