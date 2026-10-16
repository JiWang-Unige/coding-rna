#!/usr/bin/env bash
#SBATCH --job-name=M28FULLBOUNDCHK
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
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts" CUDA_VISIBLE_DEVICES=""
out=outputs/M28-FULL-GRAMMAR-BOUNDARY-CHECKS
test ! -e "$out"; mkdir "$out"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$out/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$out/STATUS"; printf '%s\n' "$SLURM_JOB_ID" > "$out/JOBID"
python - <<'PY'
import json
import os
import time
from pathlib import Path
import unittest
import torch
out = Path('outputs/M28-FULL-GRAMMAR-BOUNDARY-CHECKS')
sources = (('M28-FULL-GRAMMAR-OBSERVED-CHECKS','13427476'),
           ('M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS','13427817'))
for name, job in sources:
    root = Path('outputs')/name
    summary = json.loads((root/'summary.json').read_text())
    if (root/'STATUS').read_text().strip() != 'COMPLETED' or (root/'JOBID').read_text().strip() != job:
        raise SystemExit('Matching closed source required: '+name)
    if not summary['completed'] or summary['job_id'] != job or summary['optimizer_steps_this_run'] != 0:
        raise SystemExit('Source identity or zero-update scope differs: '+name)
torch.set_num_threads(2)
started = time.perf_counter()
with (out/'started.json').open('x') as handle:
    json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], source_jobs=[job for _,job in sources],
        synthetic_only=True, boundary_prior_chosen=False, optimizer_steps_this_run=0),handle,indent=2)
suite = unittest.defaultTestLoader.discover('tests', pattern='test_m28_full_grammar_boundary_messages.py')
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.testsRun != 10 or not result.wasSuccessful(): raise SystemExit(1)
summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
    protocol='reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md',
    source_jobs=[job for _,job in sources], implementation='native_Torch_complete_gene_boundary_message_composition',
    tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
    seconds=time.perf_counter()-started, synthetic_only=True, CPU_first_order_verified=True,
    complete_global_grammar_logZ_and_all_gradients_match_independent_and_original_oracles=True,
    every_interior_cut_and_single_base_chunks_preserve_mass_and_gradients=True,
    actual_carry_motif_stages_exon_cap_and_role_precedence_across_cuts_verified=True,
    core_observation_NLL_and_full_context_gradients_verified=True,
    incomplete_global_genes_not_accepted=True, unresolved_and_outside_genes_marginalized=True,
    boundary_prior_chosen=False, context_extent_or_padding_rule_selected=False,
    loss_normalizer_or_fit_coefficients_selected=False,
    optimizer_steps_this_run=0, source_metadata_summaries_read=True,
    research_DNA_features_labels_or_weights_opened=False, model_or_backbone_forward=False,
    original_R4_policy_changed=False, observation_rows_deleted=False,
    DEV_used=False, independent_species_used=False, A33_outputs_opened=False,
    training_or_inference_caller_changed=False,
    interpretation='artificial_exact_boundary_message_mathematics_not_actual_context_availability_resource_fit_or_performance')
with (out/'summary.json').open('x') as handle:
    json.dump(summary,handle,indent=2,allow_nan=False)
PY
printf 'COMPLETED\n' > "$out/STATUS"
