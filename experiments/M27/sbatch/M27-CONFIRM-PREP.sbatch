#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --output=outputs/M27-NATIVE-CONFIRM-R1/logs/%x_%j.out
#SBATCH --error=outputs/M27-NATIVE-CONFIRM-R1/logs/%x_%j.err
#SBATCH --job-name=M27CONFPREP
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
cd /home/users/j/jwang/coding-rna
M27_OUT="$(pwd -P)/outputs/M27-NATIVE-CONFIRM-R1"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$M27_OUT/PREP_STATUS"; fi' EXIT
printf 'RUNNING\n' > "$M27_OUT/PREP_STATUS"
python -m pytest -q -p no:cacheprovider tests/test_m27_wt_replay.py tests/test_m27_native_trace.py tests/test_m27_paired_response.py tests/test_m27_confirmation.py --junitxml="$M27_OUT/tests.xml"
printf 'COMPLETED\n' > "$M27_OUT/PREP_STATUS"
