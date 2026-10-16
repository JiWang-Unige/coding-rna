#!/usr/bin/env bash
#SBATCH --job-name=M27PREP
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=outputs/M27-ALLELIC-PREP-R1/logs/%x_%j.out
#SBATCH --error=outputs/M27-ALLELIC-PREP-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=""
OUT=outputs/M27-ALLELIC-PREP-R1
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$OUT/STATUS"
python -m pytest -q -p no:cacheprovider tests/test_m27_prepare.py --junitxml="$OUT/tests.xml"
python -u scripts/experiments/M27-ALLELIC-INTEGRITY/prepare_reference.py --data-dir data/m27_grch38_chr22 --output-dir "$OUT"
apptainer exec --cleanenv --env PYTHONDONTWRITEBYTECODE=1 --bind "$(pwd -P):/work" --pwd /work \
    refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif \
    python scripts/experiments/M27-ALLELIC-INTEGRITY/runtime_probe.py "$OUT/tiberius_runtime.json"
printf 'COMPLETED\n' > "$OUT/STATUS"
