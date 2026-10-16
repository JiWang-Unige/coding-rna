#!/usr/bin/env python3
"""New checkpoint automaton on fixed artificial 49kb cases; no research inputs."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_contig_emission_lattice_synthetic as F
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.full_grammar_viterbi import decode_full_grammar_viterbi, _forward, _trace


def main(case, out):
    torch.set_num_threads(2)
    began = time.perf_counter()
    generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (F.LENGTH,), generator=generator).tolist())
    potential_generator = torch.Generator().manual_seed(41)
    reference_root = Path('outputs/M28-A42-CONTIG-LATTICE-SCALE') / f'case_{case}'
    reference = json.loads((reference_root / 'summary.json').read_text())
    reference_job = f'13422815_{case}'
    identity = dict(case_index=case, case=F.CASES[case], length=F.LENGTH,
        windows=[list(w) for w in F.WINDOWS], strand='+', DNA_seed=0,
        potential_seed=41 if case == 1 else None,
        potential_definition='all_zero' if case == 0 else 'independent_FP32_integers_minus3_to3_all15_columns',
        context_policy='nearest_containing_genomic_window_center_per_base_lower_start_ties',
        max_exons=128, threshold=0, null=0, role_weights=1,
        motif_counts=F.S.motif_counts(sequence))
    if (not reference['completed'] or reference['job_id'] != reference_job
            or (reference_root / 'JOBID').read_text().strip() != reference_job
            or (reference_root / 'STATUS').read_text().strip() != 'COMPLETED'
            or any(reference[key] != value for key, value in identity.items())):
        raise ValueError('Fixed completed independent Python reference or artificial input identity changed')
    job = os.environ['SLURM_ARRAY_JOB_ID'] + '_' + os.environ['SLURM_ARRAY_TASK_ID']
    started = dict(**identity, job_id=job, started_utc=F.S.utc_now(), synthetic_only=True,
        backend='checkpointed_complete_boundary_state_automaton', checkpoint_bp=8192,
        score_accumulation='FP64_per_base_forward_order',
        tie_policy='background_continuation_then_state_table_arrival',
        independent_reference_job_id=reference_job,
        comparison='optimal_integer_score_only_not_tied_chain_identity_or_old_counters')
    with (out / 'started.json').open('x') as handle:
        json.dump(started, handle, indent=2, allow_nan=False)
    print(json.dumps(started, allow_nan=False), flush=True)

    def records():
        for i, (a, b) in enumerate(F.WINDOWS):
            logits = (torch.zeros(b-a, 15) if case == 0 else
                torch.randint(-3, 4, (b-a, 15), generator=potential_generator).float())
            yield dict(window_index=i, sequence=sequence[a:b], local_potentials=logits)

    joined = stitch_contig_emissions(sequence, F.WINDOWS, '+', records())
    logits = joined['local_potentials']
    upper = float((logits.double() - logits[:, :1].double()).max(dim=1).values.sum())
    segments = [list(v) for v in joined['genomic_context_segments']]
    if upper != reference['run']['pointwise_upper_bound'] or segments != reference['run']['genomic_context_segments']:
        raise ValueError('Fixed integer context preparation differs from completed reference')
    preparation_seconds = time.perf_counter() - began
    print(json.dumps(dict(stage='new_automaton_decode_started', utc=F.S.utc_now(),
        pointwise_upper_bound=upper, preparation_seconds=preparation_seconds)), flush=True)
    began = time.perf_counter()
    result = decode_full_grammar_viterbi(sequence, logits, 128)
    seconds = time.perf_counter() - began
    F.S.verify_selected(sequence, logits, result, 128)
    if result['gain'] != reference['run']['gain']:
        raise ValueError('New complete-set max score differs from independent completed Python optimum')
    if not 0 <= result['gain'] <= upper or (case == 0 and (result['gain'] != 0 or result['chains'])):
        raise ValueError('Pointwise score bound or background-first zero null failed')
    counts = result['counts']; states = counts['state_count']
    if (counts['checkpoint_bytes'] != (F.LENGTH // 8192 + 1) * states * 8
            or counts['maximum_block_predecessor_bytes'] != 8192 * states * 4
            or counts['input_logit_bytes'] != F.LENGTH * 15 * 4
            or counts['forward_and_recomputed_bases'] != 2 * F.LENGTH):
        raise ValueError('Measured checkpoint representation differs from its explicit arrays')
    if not _forward.nopython_signatures or not _trace.nopython_signatures:
        raise ValueError('Both new kernels must actually execute in nopython mode')
    run = dict(decode_seconds_including_cold_JIT=seconds, preparation_seconds=preparation_seconds,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        gain=result['gain'], independent_reference_gain=reference['run']['gain'],
        pointwise_upper_bound=upper, selected_count=len(result['chains']),
        selected_chains=result['chains'], chain_gains=result['chain_gains'], counts=counts,
        genomic_context_segments=segments,
        forward_nopython_signatures=len(_forward.nopython_signatures),
        trace_nopython_signatures=len(_trace.nopython_signatures),
        selected_DNA_ORF_motifs_carry_cap_nonoverlap_and_integer_scores_recomputed=True,
        optimal_integer_score_matches_independent_Python_reference=True)
    with (out / 'records.jsonl').open('x') as handle:
        handle.write(json.dumps(run, allow_nan=False) + '\n')
    summary = dict(**started, completed=True, completed_utc=F.S.utc_now(), run=run,
        research_inputs_opened=False, model_or_backbone_forward=False, optimizer_steps=0,
        labels_or_catalog_opened=False, DEV_used=False, independent_species_used=False,
        A33_outputs_opened=False, stopped_A43_endpoint_prototype_used=False,
        old_independent_Python_solver_rerun=False, A44_execution_chain_reopened=False,
        interpretation='artificial_49kb_resource_and_independent_optimum_not_chromosome_or_annotation_performance')
    with (out / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps({k:v for k,v in run.items() if k not in ('selected_chains', 'chain_gains')}, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=int, choices=(0,1), required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    main(args.case, args.out)
