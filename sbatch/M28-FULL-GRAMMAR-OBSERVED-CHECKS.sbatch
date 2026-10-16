#!/usr/bin/env bash
#SBATCH --job-name=M28FULLOBSCHK
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
export PYTHONDONTWRITEBYTECODE=1 NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD:$PWD/scripts/experiments/M28-METHOD-RESTART" CUDA_VISIBLE_DEVICES=""
out=outputs/M28-FULL-GRAMMAR-OBSERVED-CHECKS
test ! -e "$out"; mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python - <<'PY'
import json
import os
from pathlib import Path
import unittest
import torch
from src.m28.full_grammar_partition_compiled_reverse import _accumulate_reverse_nodes
torch.set_num_threads(2)
suite = unittest.defaultTestLoader.discover('tests', pattern='test_m28_full_grammar_observed_partition.py')
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.testsRun != 10 or not result.wasSuccessful(): raise SystemExit(1)
if not _accumulate_reverse_nodes.nopython_signatures: raise SystemExit('No actual nopython reverse')
summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
    protocol='full_grammar_partition_preparation.md',
    original_compiled_reverse_checks_job='13427223', original_complete_window_job='13427267',
    implementation='partial_role_event_complete_gene_sets_with_compiled_reverse',
    tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
    synthetic_only=True, CPU_first_order_only=True,
    actual_nopython_signature_count=len(_accumulate_reverse_nodes.nopython_signatures),
    artificial_observed_mass_and_NLL_full_gradients_verified=True,
    original_unconditional_forward_and_counts_on_all_unknown_verified=True,
    fine_coarse_unknown_padding_mapping_verified_on_artificial_targets=True,
    background_and_intron_observation_constraints_verified=True,
    infeasible_observation_raises_without_boundary_repair=True,
    observed_likelihood_implemented=True, loss_normalizer_or_fit_coefficients_selected=False,
    optimizer_steps_this_run=0, actual_research_inputs_opened=False,
    model_or_backbone_forward=False, research_candidate_search=False,
    original_R4_targets_or_catalog_opened=False, DEV_used=False, independent_species_used=False,
    training_or_inference_caller_changed=False,
    interpretation='artificial_partial_observation_likelihood_not_actual_R4_feasibility_fit_or_performance')
with Path('outputs/M28-FULL-GRAMMAR-OBSERVED-CHECKS/summary.json').open('x') as handle:
    json.dump(summary, handle, indent=2, allow_nan=False)
PY
printf 'COMPLETED\n' > "$out/STATUS"
