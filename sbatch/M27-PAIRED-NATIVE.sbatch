#!/usr/bin/env bash
#SBATCH --job-name=M27PAIRED
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=outputs/M27-PAIRED-NATIVE-R1/logs/%x_%j.out
#SBATCH --error=outputs/M27-PAIRED-NATIVE-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
M27_OUT="$(pwd -P)/outputs/M27-PAIRED-NATIVE-R1"
test "$(cat "$M27_OUT/PREP_STATUS")" = COMPLETED
test ! -e "$M27_OUT/STATUS"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$M27_OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$M27_OUT/STATUS"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate annevo
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=8
python -c 'import torch; print("Torch",torch.__version__,"CUDA",torch.version.cuda,"devices",torch.cuda.device_count()); assert torch.cuda.is_available() and torch.cuda.device_count()==1'
python -u scripts/experiments/M27-ALLELIC-INTEGRITY/paired_run.py --output-dir "$M27_OUT"
printf 'COMPLETED\n' > "$M27_OUT/STATUS"
