#!/bin/bash
#SBATCH --job-name=M24-DIRECT-STRUCTURE-DIAGNOSTIC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC/logs/%x_%j.out
#SBATCH --error=outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC/logs/%x_%j.err

set -euo pipefail

cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate coding-rna
set -u

python scripts/eval_structure_diagnostic.py \
  --repo-root . \
  --out-dir outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC

touch outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC/COMPLETED
