#!/usr/bin/env bash
# Create the dedicated `annevo` conda env (sanctioned exception to coding-rna: ANNEVO needs
# Py3.10/torch2.1/cu12.1, incompatible with the project env). Run detached on the baobab
# login node (needs internet for conda channels). Poll outputs/SETUP-ANNEVO-ENV/setup.log.
set -uo pipefail
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
cd /srv/beegfs/scratch/shares/ds4dh/common/coding-rna

echo "START=$(date +%s) $(date)"
echo "=== mamba env create -n annevo from ANNEVO.yml ==="
if mamba env create -f refs/repos/annevo-2026/ANNEVO.yml -n annevo; then
  echo "CREATE_EXIT=0 (mamba)"
else
  echo "MAMBA_FAILED rc=$? — retrying with conda"
  conda env create -f refs/repos/annevo-2026/ANNEVO.yml -n annevo
  echo "CREATE_EXIT=$? (conda)"
fi

echo "=== sanity: torch + cuda visibility (CPU import only on login node) ==="
conda activate annevo && python -c "import torch, h5py, numba; print('torch', torch.__version__, 'cuda_build', torch.version.cuda)" && echo "IMPORT_OK"
echo "END=$(date +%s) $(date)"
echo "SETUP_DONE"
