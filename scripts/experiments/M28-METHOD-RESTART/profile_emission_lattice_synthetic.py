#!/usr/bin/env python3
"""Fixed synthetic full-window A42 resource check, no research inputs."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import time

import torch

from src.m28.emission_lattice import decode_emission_lattice
from src.m28.phase_evidence_chain import candidate_roles


WIDTH = 24576
STOP = {'TAA', 'TAG', 'TGA'}
CASES = ('known_128_exons', 'iid_uniform_zero_potentials')


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def known_fixture():
    sequence, chain = 'C' * 512, []
    for index in range(128):
        exon = 'ATGAAA' if index == 0 else 'CCCTAA' if index == 127 else 'CCC'
        a = len(sequence)
        sequence += exon
        chain.append((a, len(sequence)))
        if index != 127:
            sequence += 'GTTAG'
    sequence += 'C' * (WIDTH - len(sequence))
    chain = tuple(chain)
    logits = torch.full((WIDTH, 15), -4., dtype=torch.float64)
    logits[:, 0] = 0.
    roles = candidate_roles(chain, WIDTH)
    positions = torch.arange(chain[0][0], chain[-1][1])
    logits[positions, roles[positions]] = 2.
    # Every base in the target span has maximum role margin 2, elsewhere 0.
    upper_bound = 2. * (chain[-1][1] - chain[0][0])
    return sequence, logits, chain, upper_bound


def random_fixture():
    generator = torch.Generator().manual_seed(0)
    alphabet = 'ACGT'
    sequence = ''.join(alphabet[i] for i in torch.randint(0, 4, (WIDTH,), generator=generator).tolist())
    return sequence, torch.zeros(WIDTH, 15, dtype=torch.float64)


def verify_selected(sequence, logits, result, cap):
    """Independent concatenated-CDS legality and direct role arithmetic."""
    if result['reference_used'] is not False or len(result['chains']) != len(result['chain_gains']):
        raise ValueError('Reference-free complete chain/gain pairs required')
    previous_gene_end = 0
    for chain, gain in zip(result['chains'], result['chain_gains']):
        if not 0 < len(chain) <= cap or chain[0][0] < previous_gene_end:
            raise ValueError('Complete capped, nonoverlapping gene spans required')
        if chain[0][1] - chain[0][0] < 3 or chain[-1][1] - chain[-1][0] < 3:
            raise ValueError('Contiguous terminal codons required')
        coding = ''
        previous_exon_end = None
        for a, b in chain:
            if not 0 <= a < b <= len(sequence):
                raise ValueError('Valid oriented exon coordinates required')
            if previous_exon_end is not None:
                if (a - previous_exon_end < 5 or sequence[previous_exon_end:previous_exon_end + 2] not in ('GT', 'GC')
                        or sequence[a - 2:a] != 'AG'):
                    raise ValueError('Original splice motif and minimum intron grammar required')
            coding += sequence[a:b]
            previous_exon_end = b
        codons = [coding[i:i + 3] for i in range(0, len(coding) - 2, 3)]
        if (any(base not in 'ACGT' for base in coding) or len(coding) % 3 or
                not coding.startswith('ATG') or coding[-3:] not in STOP or
                any(codon in STOP for codon in codons[:-1])):
            raise ValueError('Actual cross-exon ORF grammar differs')
        roles = candidate_roles(chain, len(sequence))
        positions = torch.arange(chain[0][0], chain[-1][1])
        direct = float((logits[positions, roles[positions]].double() - logits[positions, 0].double()).sum())
        if direct != gain or gain <= 0.:
            raise ValueError('Synthetic integer FP64 role score differs or is nonpositive')
        previous_gene_end = chain[-1][1]
    if sum(result['chain_gains']) != result['gain']:
        raise ValueError('Selected complete chain gain sum differs')


def motif_counts(sequence):
    return dict(start=sum(sequence[i:i + 3] == 'ATG' for i in range(len(sequence) - 2)),
                donor=sum(sequence[i:i + 2] in ('GT', 'GC') for i in range(1, len(sequence) - 1)),
                acceptor=sum(sequence[i:i + 2] == 'AG' for i in range(len(sequence) - 2)),
                stop=sum(sequence[i:i + 3] in STOP for i in range(len(sequence) - 2)))


def main(case, out):
    torch.set_num_threads(2)
    expected = None
    if case == 0:
        sequence, logits, expected, upper_bound = known_fixture()
        caps = (128, 127)
    else:
        sequence, logits = random_fixture()
        upper_bound, caps = 0., (128,)
    if len(sequence) != WIDTH:
        raise ValueError('Actual full window length required')
    started = dict(case=CASES[case],case_index=case,width=WIDTH,exon_caps=list(caps),
        synthetic_only=True,seed=0 if case == 1 else None,
        potential_definition='all_zero' if case == 1 else 'target_roles_2_other_roles_minus4_intergenic0',
        motif_counts=motif_counts(sequence),pointwise_upper_bound=upper_bound,started_utc=utc_now(),
        job_id=os.environ['SLURM_ARRAY_JOB_ID']+'_'+os.environ['SLURM_ARRAY_TASK_ID'])
    with (out/'started.json').open('x') as handle:
        json.dump(started, handle, ensure_ascii=False, allow_nan=False, indent=2)
    print(json.dumps(started, allow_nan=False), flush=True)
    records = []
    with (out/'records.jsonl').open('x') as handle:
        for cap in caps:
            print(json.dumps(dict(stage='decode_started',exon_cap=cap,utc=utc_now())), flush=True)
            began = time.perf_counter()
            result = decode_emission_lattice(sequence, logits, cap)
            seconds = time.perf_counter() - began
            verify_selected(sequence, logits, result, cap)
            if case == 0 and cap == 128:
                if result['chains'] != (expected,) or result['gain'] != upper_bound:
                    raise ValueError('Known 128-exon optimum must attain the independent pointwise upper bound')
            elif case == 0:
                if expected in result['chains'] or result['gain'] >= upper_bound:
                    raise ValueError('Exon cap must exclude the unique 128-exon pointwise optimum')
            elif result['chains'] or result['gain'] != 0.:
                raise ValueError('Zero potentials must return null regardless of grammar density')
            record = dict(exon_cap=cap,seconds=seconds,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                selected_chains=result['chains'],chain_gains=result['chain_gains'],gain=result['gain'],
                counts=result['counts'],actual_selected_ORF_roles_and_span_recomputed=True,
                known_upper_bound_and_cap_or_zero_null_verified=True)
            handle.write(json.dumps(record,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps({k:v for k,v in record.items() if k not in ('selected_chains','chain_gains')},allow_nan=False),flush=True)
            records.append(record)
    summary = dict(**started,completed=True,completed_utc=utc_now(),runs=records,
        optimizer_steps=0,research_inputs_opened=False,model_or_backbone_forward=False,
        new_research_candidate_search=False,DEV_used=False,A33_outputs_opened=False,
        interpretation='synthetic_resource_and_known_answer_check_not_annotation_performance_or_genome_speed')
    with (out/'summary.json').open('x') as handle:
        json.dump(summary, handle, ensure_ascii=False, allow_nan=False, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=int,choices=(0,1),required=True)
    parser.add_argument('--out',type=Path,required=True)
    arguments = parser.parse_args()
    main(arguments.case,arguments.out)
