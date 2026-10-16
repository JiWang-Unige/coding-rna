#!/usr/bin/env bash
# Install the DL/data stack INTO the project env `coding-rna` (Py3.12) for our own from-scratch
# screen-anchor training harness. Project-env discipline: our model code lives in coding-rna
# (the `annevo` env is the ANNEVO baseline tool only). Run detached on the baobab login node
# (needs internet). Poll outputs/SETUP-CODINGRNA-TORCH/setup.log.
set -uo pipefail
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
cd /srv/beegfs/scratch/shares/ds4dh/common/coding-rna

echo "START=$(date +%s) $(date)"
echo "=== mamba install torch(cuda12.1)+numpy into coding-rna ==="
# FASTA is parsed in pure python in the harness (no pyfaidx dep, which lives on bioconda).
mamba install -n coding-rna -y -c pytorch -c nvidia -c conda-forge \
  pytorch pytorch-cuda=12.1 numpy
echo "INSTALL_EXIT=$?"

echo "=== sanity import ==="
set +u; conda activate coding-rna; set -u
python -c "import torch, numpy; print('torch', torch.__version__, 'cuda_build', torch.version.cuda, '| numpy', numpy.__version__)" && echo "IMPORT_OK"
echo "END=$(date +%s) $(date)"
echo "SETUP_DONE"
