#!/usr/bin/env python3
"""Bounded artificial three-context joint decode scale probe; no research data."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import torch

import profile_emission_lattice_synthetic as S
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice import decode_emission_lattice


LENGTH = 49152
WINDOWS = ((0, 24576), (12288, 36864), (24576, 49152))
CASES = ('iid_uniform_zero_potentials', 'iid_uniform_integer_window_potentials')


def main(case, out):
    torch.set_num_threads(2)
    began = time.perf_counter()
    dna_generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (LENGTH,), generator=dna_generator).tolist())
    potential_generator = torch.Generator().manual_seed(41)
    job_id = os.environ['SLURM_ARRAY_JOB_ID'] + '_' + os.environ['SLURM_ARRAY_TASK_ID']
    started = dict(case_index=case,case=CASES[case],length=LENGTH,windows=WINDOWS,strand='+',
        DNA_seed=0,potential_seed=41 if case == 1 else None,
        potential_definition='all_zero' if case == 0 else 'independent_FP32_integers_minus3_to3_all15_columns',
        context_policy='nearest_containing_genomic_window_center_per_base_lower_start_ties',
        max_exons=128,threshold=0,null=0,role_weights=1,score_accumulation='FP64',
        motif_counts=S.motif_counts(sequence),synthetic_only=True,job_id=job_id,started_utc=S.utc_now())
    with (out / 'started.json').open('x') as handle:
        json.dump(started,handle,indent=2,allow_nan=False)
    print(json.dumps(started,allow_nan=False),flush=True)

    def records():
        for i, (a, b) in enumerate(WINDOWS):
            logits = (torch.zeros(b - a, 15) if case == 0 else
                      torch.randint(-3, 4, (b - a, 15), generator=potential_generator).float())
            yield dict(window_index=i,sequence=sequence[a:b],local_potentials=logits)

    joined = stitch_contig_emissions(sequence,WINDOWS,'+',records())
    preparation_seconds = time.perf_counter() - began
    logits = joined['local_potentials']
    upper_bound = float((logits.double() - logits[:, :1].double()).max(dim=1).values.sum())
    print(json.dumps(dict(stage='joint_decode_started',utc=S.utc_now(),
        preparation_seconds=preparation_seconds,pointwise_upper_bound=upper_bound,
        genomic_context_segments=joined['genomic_context_segments'])),flush=True)
    began = time.perf_counter()
    result = decode_emission_lattice(joined['sequence'],logits,128)
    decode_seconds = time.perf_counter() - began
    S.verify_selected(joined['sequence'],logits,result,128)
    if not 0 <= result['gain'] <= upper_bound:
        raise ValueError('Nonoverlapping complete gene set violates the pointwise role upper bound')
    if case == 0 and (result['chains'] or result['gain'] != 0):
        raise ValueError('Dense all-zero input must return null without a zero shortcut')
    record = dict(decode_seconds=decode_seconds,preparation_seconds=preparation_seconds,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        pointwise_upper_bound=upper_bound,selected_chains=result['chains'],
        chain_gains=result['chain_gains'],selected_count=len(result['chains']),gain=result['gain'],
        counts=result['counts'],genomic_context_segments=joined['genomic_context_segments'],
        actual_selected_DNA_ORF_motifs_roles_cap_and_spans_recomputed=True,
        actual_direct_integer_FP64_chain_and_total_scores_recomputed=True,
        pointwise_upper_bound_and_zero_null_when_applicable_verified=True)
    with (out / 'records.jsonl').open('x') as handle:
        handle.write(json.dumps(record,allow_nan=False) + '\n')
    summary = dict(**started,completed=True,completed_utc=S.utc_now(),run=record,
        optimizer_steps=0,research_inputs_opened=False,model_or_backbone_forward=False,
        new_research_candidate_search=False,labels_or_catalog_opened=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,
        interpretation='artificial_resource_probe_not_annotation_performance_or_chromosome_feasibility')
    with (out / 'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in record.items() if k not in ('selected_chains','chain_gains')},allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=int,choices=(0,1),required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    main(args.case,args.out)
