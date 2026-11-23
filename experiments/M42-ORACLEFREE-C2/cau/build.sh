#!/bin/bash
# CPU job: build causal-conv1d 1.6.0 and mamba-ssm 2.3.0 from source against the pip torch 2.5.1 (cu124) of $G/env, with the env's CUDA 12.4 nvcc and gcc 12
# (the prebuilt GitHub wheels used on ont need GLIBC 2.32; this cluster has 2.28).  Same package versions as ont's plantcad2 env.
umask 077; set -uo pipefail
G=$HOME/.cache/fontconfig-x2/g10; E=$G/env; mkdir -p $G/tmp
export PATH=$E/bin:$PATH CUDA_HOME=$E TMPDIR=$G/tmp PIP_CACHE_DIR=$G/pipcache MAX_JOBS=${SLURM_CPUS_PER_TASK:-8} \
       MAMBA_FORCE_BUILD=TRUE CAUSAL_CONV1D_FORCE_BUILD=TRUE MAMBA_SKIP_CUDA_BUILD=FALSE CAUSAL_CONV1D_SKIP_CUDA_BUILD=FALSE \
       CC=$E/bin/x86_64-conda-linux-gnu-gcc CXX=$E/bin/x86_64-conda-linux-gnu-g++ CUDAHOSTCXX=$E/bin/x86_64-conda-linux-gnu-g++
echo "== $(date)"; nvcc --version | tail -2; $CXX --version | head -1
cd $G/src
for p in causal_conv1d-1.6.0.tar.gz mamba_ssm-2.3.0.tar.gz; do
  echo "== build $p $(date)"
  $E/bin/python -m pip install -v --no-build-isolation --no-deps --no-index $p || { echo "BUILD_FAIL $p"; exit 1; }
done
$E/bin/python - <<'P'
import importlib, os, torch
for m in ('causal_conv1d_cuda', 'selective_scan_cuda'):
    x = importlib.import_module(m); print(m, os.path.getsize(x.__file__))
import mamba_ssm, causal_conv1d; print('mamba_ssm', mamba_ssm.__version__, 'causal_conv1d', causal_conv1d.__version__)
P
echo "BUILD_DONE $(date)"
