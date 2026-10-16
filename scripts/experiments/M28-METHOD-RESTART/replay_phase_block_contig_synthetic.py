#!/usr/bin/env python3
"""Replay CLOSED artificial cases with same-target phase-block array dispatch."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_contig_emission_lattice_synthetic as P
import profile_emission_lattice_synthetic as S
from replay_batched_contig_synthetic import COUNT_PAIRS
from replay_indexed_contig_synthetic import ANSWER_FIELDS, IDENTITY_FIELDS, SOURCE, SOURCE_ARRAY
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice_phase_block import decode_phase_block_emission_lattice


CHECKS = Path('outputs/M28-A42-PHASE-BLOCK-CHECKS')
CHECKS_JOB = '13423246'
PREVIOUS = Path('outputs/M28-A42-BATCHED-CONTIG-REPLAY')
PREVIOUS_ARRAY = '13423130'
SCOPE_FIELDS = ('research_inputs_opened', 'model_or_backbone_forward', 'new_research_candidate_search',
                'labels_or_catalog_opened', 'DEV_used', 'independent_species_used', 'A33_outputs_opened')


def main(case, out):
    torch.set_num_threads(2)
    source, source_job = SOURCE / ('case_' + str(case)), SOURCE_ARRAY + '_' + str(case)
    previous, previous_job = PREVIOUS / ('case_' + str(case)), PREVIOUS_ARRAY + '_' + str(case)
    for root, expected in ((CHECKS, CHECKS_JOB), (source, source_job), (previous, previous_job)):
        if root.joinpath('STATUS').read_text().strip() != 'COMPLETED' or root.joinpath('JOBID').read_text().strip() != expected:
            raise ValueError('Assigned checks and original artificial cases must be closed first')
    dna_generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (P.LENGTH,), generator=dna_generator).tolist())
    potential_generator = torch.Generator().manual_seed(41)
    job_id = os.environ['SLURM_ARRAY_JOB_ID'] + '_' + os.environ['SLURM_ARRAY_TASK_ID']
    started = dict(case_index=case, case=P.CASES[case], length=P.LENGTH, windows=P.WINDOWS, strand='+',
        DNA_seed=0, potential_seed=41 if case == 1 else None,
        potential_definition='all_zero' if case == 0 else 'independent_FP32_integers_minus3_to3_all15_columns',
        context_policy='nearest_containing_genomic_window_center_per_base_lower_start_ties',
        max_exons=128, threshold=0, null=0, role_weights=1, score_accumulation='FP64',
        motif_counts=S.motif_counts(sequence), synthetic_only=True, job_id=job_id,
        source_job=source_job, previous_batch_job=previous_job, checks_job=CHECKS_JOB, started_utc=S.utc_now())
    with (out / 'started.json').open('x') as handle:
        json.dump(started, handle, indent=2, allow_nan=False)

    def records():
        for i, (a, b) in enumerate(P.WINDOWS):
            logits = (torch.zeros(b - a, 15) if case == 0 else
                      torch.randint(-3, 4, (b - a, 15), generator=potential_generator).float())
            yield dict(window_index=i, sequence=sequence[a:b], local_potentials=logits)

    joined = stitch_contig_emissions(sequence, P.WINDOWS, '+', records())
    logits = joined['local_potentials']
    upper_bound = float((logits.double() - logits[:, :1].double()).max(dim=1).values.sum())
    print(json.dumps(dict(**started, stage='phase_block_decode_started', pointwise_upper_bound=upper_bound), allow_nan=False), flush=True)
    began = time.perf_counter()
    result = decode_phase_block_emission_lattice(joined['sequence'], logits, 128)
    seconds = time.perf_counter() - began
    S.verify_selected(joined['sequence'], logits, result, 128)
    if not 0 <= result['gain'] <= upper_bound or (case == 0 and (result['chains'] or result['gain'] != 0)):
        raise ValueError('Original upper bound or zero-null contract differs')
    counts = result['counts']
    if not (counts['state_update_proposals'] <= counts['logical_exon_extensions'] and
            counts['long_array_blocks'] <= counts['long_nonempty_history_groups'] <= counts['long_history_groups']):
        raise ValueError('Logical/physical and nonempty/block counters differ')
    record = dict(decode_seconds=seconds, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        pointwise_upper_bound=upper_bound, selected_chains=result['chains'], chain_gains=result['chain_gains'],
        selected_count=len(result['chains']), gain=result['gain'], counts=counts,
        genomic_context_segments=joined['genomic_context_segments'],
        actual_selected_DNA_ORF_motifs_roles_cap_and_spans_recomputed=True,
        actual_direct_integer_FP64_chain_and_total_scores_recomputed=True,
        pointwise_upper_bound_and_zero_null_when_applicable_verified=True)
    with (out / 'natural_records.jsonl').open('x') as handle:
        handle.write(json.dumps(record, allow_nan=False) + '\n')
    # Both prior ANSWERS are opened only after the new natural output closes.
    original = json.loads(source.joinpath('summary.json').read_text())
    batch = json.loads(previous.joinpath('summary.json').read_text())
    identity, actual = json.loads(json.dumps(started)), json.loads(json.dumps(record))
    for saved, expected in ((original, source_job), (batch, previous_job)):
        if (saved['completed'] is not True or saved['job_id'] != expected or
                any(saved[key] != identity[key] for key in IDENTITY_FIELDS) or
                any(saved[key] is not False for key in SCOPE_FIELDS) or saved['optimizer_steps'] != 0):
            raise ValueError('Prior synthetic configuration or scope identity differs')
        if any(actual[key] != saved['run'][key] for key in ANSWER_FIELDS if key != 'counts'):
            raise ValueError('Phase-block chains/order or FP64 gains differ from a prior answer')
    if any(counts[new] != original['run']['counts'][old] for new, old in COUNT_PAIRS):
        raise ValueError('Original logical extensions/terminal/cap/active differs')
    if any(counts[key] != value for key, value in batch['run']['counts'].items()):
        raise ValueError('Previous batch physical or logical counts differ')
    summary = dict(**started, completed=True, completed_utc=S.utc_now(), run=record,
        actual_saved_synthetic_configuration_and_both_prior_jobs_identity_verified=True,
        actual_new_natural_search_saved_and_closed_before_both_prior_answers_opened=True,
        actual_all_selected_chain_identities_order_and_FP64_gains_equal_both_prior_answers=True,
        actual_original_logical_and_all_previous_batch_counts_equal=True,
        original_decode_seconds=original['run']['decode_seconds'], previous_batch_decode_seconds=batch['run']['decode_seconds'],
        optimizer_steps=0, research_inputs_opened=False, model_or_backbone_forward=False,
        new_research_candidate_search=False, labels_or_catalog_opened=False, DEV_used=False,
        independent_species_used=False, A33_outputs_opened=False,
        interpretation='same_artificial_case_not_fair_speed_or_annotation_performance_nonempty_groups_are_not_total_groups')
    with (out / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps({key: value for key, value in summary.items() if key != 'run'}, allow_nan=False), flush=True)
    print(json.dumps({key: value for key, value in record.items() if key not in ('selected_chains', 'chain_gains')}, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=int, choices=(0, 1), required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    main(args.case, args.out)
