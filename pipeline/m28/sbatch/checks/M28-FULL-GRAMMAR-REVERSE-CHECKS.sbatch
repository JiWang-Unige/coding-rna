#!/usr/bin/env bash
#SBATCH --job-name=M28FULLZREVCHK
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:02:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/M28-PILOT-R5/logs/%x_%j.out
#SBATCH --error=outputs/M28-PILOT-R5/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate generanno; set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts" CUDA_VISIBLE_DEVICES=""
out=outputs/M28-FULL-GRAMMAR-REVERSE-CHECKS
test ! -e "$out"; mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python - <<'PY'
import json
import os
from pathlib import Path
import unittest
import torch
torch.set_num_threads(2)
suite = unittest.defaultTestLoader.discover('tests', pattern='test_m28_full_grammar_partition_reverse.py')
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.testsRun != 9 or not result.wasSuccessful(): raise SystemExit(1)
summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
    protocol='full_grammar_partition_preparation.md', original_mathematical_checks_job='13426777',
    indexed_mathematical_checks_job='13427003',
    implementation='compact_grammar_specific_explicit_first_order_CPU_reverse', tests_run=result.testsRun,
    failures=len(result.failures), errors=len(result.errors), synthetic_only=True, higher_order_gradients_supported=False, GPU_supported=False,
    mathematical_partition_and_gradient_checks_passed=True, optimizer_steps_this_run=0,
    actual_research_inputs_opened=False, model_or_backbone_forward=False,
    research_candidate_search=False, labels_or_catalog_opened=False,
    DEV_used=False, independent_species_used=False,
    observed_mask_or_gold_likelihood_implemented=False,
    training_or_inference_caller_changed=False,
    interpretation='artificial_full_grammar_logZ_and_gradients_not_training_or_annotation_performance')
with Path('outputs/M28-FULL-GRAMMAR-REVERSE-CHECKS/summary.json').open('x') as handle:
    json.dump(summary, handle, indent=2, allow_nan=False)
PY
printf 'COMPLETED\n' > "$out/STATUS"
