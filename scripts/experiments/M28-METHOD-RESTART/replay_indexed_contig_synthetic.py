#!/usr/bin/env python3
"""Replay two CLOSED artificial cases with the isolated indexed A42 decoder."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_contig_emission_lattice_synthetic as P
import profile_emission_lattice_synthetic as S
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice_indexed import decode_indexed_emission_lattice


SOURCE = Path('outputs/M28-A42-CONTIG-LATTICE-SCALE')
SOURCE_ARRAY = '13422815'
IDENTITY_FIELDS = ('case_index','case','length','windows','strand','DNA_seed','potential_seed',
                   'potential_definition','context_policy','max_exons','threshold','null',
                   'role_weights','score_accumulation','motif_counts','synthetic_only')
ANSWER_FIELDS = ('selected_chains','chain_gains','selected_count','gain','counts',
                 'pointwise_upper_bound','genomic_context_segments')


def main(case, out):
    torch.set_num_threads(2)
    source = SOURCE / ('case_' + str(case))
    source_job = SOURCE_ARRAY + '_' + str(case)
    if (source.joinpath('STATUS').read_text().strip() != 'COMPLETED' or
            source.joinpath('JOBID').read_text().strip() != source_job):
        raise ValueError('The assigned original artificial case must be closed with its original job ID')
    dna_generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (P.LENGTH,), generator=dna_generator).tolist())
    potential_generator = torch.Generator().manual_seed(41)
    job_id = os.environ['SLURM_ARRAY_JOB_ID'] + '_' + os.environ['SLURM_ARRAY_TASK_ID']
    started = dict(case_index=case,case=P.CASES[case],length=P.LENGTH,windows=P.WINDOWS,strand='+',
        DNA_seed=0,potential_seed=41 if case == 1 else None,
        potential_definition='all_zero' if case == 0 else 'independent_FP32_integers_minus3_to3_all15_columns',
        context_policy='nearest_containing_genomic_window_center_per_base_lower_start_ties',
        max_exons=128,threshold=0,null=0,role_weights=1,score_accumulation='FP64',
        motif_counts=S.motif_counts(sequence),synthetic_only=True,job_id=job_id,
        source_job=source_job,started_utc=S.utc_now())
    with (out / 'started.json').open('x') as handle:
        json.dump(started,handle,indent=2,allow_nan=False)

    def records():
        for i, (a, b) in enumerate(P.WINDOWS):
            logits = (torch.zeros(b - a, 15) if case == 0 else
                      torch.randint(-3, 4, (b - a, 15), generator=potential_generator).float())
            yield dict(window_index=i,sequence=sequence[a:b],local_potentials=logits)

    joined = stitch_contig_emissions(sequence,P.WINDOWS,'+',records())
    logits = joined['local_potentials']
    upper_bound = float((logits.double() - logits[:, :1].double()).max(dim=1).values.sum())
    print(json.dumps(dict(**started,stage='indexed_decode_started',pointwise_upper_bound=upper_bound),allow_nan=False),flush=True)
    began = time.perf_counter()
    result = decode_indexed_emission_lattice(joined['sequence'],logits,128)
    seconds = time.perf_counter() - began
    S.verify_selected(joined['sequence'],logits,result,128)
    if not 0 <= result['gain'] <= upper_bound or (case == 0 and (result['chains'] or result['gain'] != 0)):
        raise ValueError('Original pointwise upper bound or zero-null contract differs')
    record = dict(decode_seconds=seconds,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        pointwise_upper_bound=upper_bound,selected_chains=result['chains'],chain_gains=result['chain_gains'],
        selected_count=len(result['chains']),gain=result['gain'],counts=result['counts'],
        genomic_context_segments=joined['genomic_context_segments'],
        actual_selected_DNA_ORF_motifs_roles_cap_and_spans_recomputed=True,
        actual_direct_integer_FP64_chain_and_total_scores_recomputed=True,
        pointwise_upper_bound_and_zero_null_when_applicable_verified=True)
    with (out / 'natural_records.jsonl').open('x') as handle:
        handle.write(json.dumps(record,allow_nan=False) + '\n')
    # Original ANSWERS are read only after the new natural search is saved/closed.
    original = json.loads(source.joinpath('summary.json').read_text())
    identity = json.loads(json.dumps(started,allow_nan=False))
    if (original['completed'] is not True or original['job_id'] != source_job or
            any(original[key] != identity[key] for key in IDENTITY_FIELDS) or
            any(original[key] is not False for key in ('research_inputs_opened','model_or_backbone_forward',
                'new_research_candidate_search','labels_or_catalog_opened','DEV_used',
                'independent_species_used','A33_outputs_opened')) or original['optimizer_steps'] != 0):
        raise ValueError('Original synthetic input/configuration or scope identity differs')
    actual = json.loads(json.dumps(record,allow_nan=False))
    if any(actual[key] != original['run'][key] for key in ANSWER_FIELDS):
        raise ValueError('Indexed decoder differs from the closed original artificial chains, gains or counters')
    summary = dict(**started,completed=True,completed_utc=S.utc_now(),run=record,
        actual_saved_synthetic_configuration_identity_verified=True,
        actual_new_natural_search_saved_and_closed_before_original_answers_opened=True,
        actual_all_selected_chains_FP64_gains_and_all_counters_equal_original=True,
        original_decode_seconds=original['run']['decode_seconds'],
        optimizer_steps=0,research_inputs_opened=False,model_or_backbone_forward=False,
        labels_or_catalog_opened=False,DEV_used=False,independent_species_used=False,A33_outputs_opened=False,
        interpretation='same_artificial_case_regression_not_fair_timing_benchmark_or_annotation_performance')
    with (out / 'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({key:value for key,value in summary.items() if key != 'run'},allow_nan=False),flush=True)
    print(json.dumps({key:value for key,value in record.items() if key not in ('selected_chains','chain_gains')},allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=int,choices=(0,1),required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    main(args.case,args.out)
