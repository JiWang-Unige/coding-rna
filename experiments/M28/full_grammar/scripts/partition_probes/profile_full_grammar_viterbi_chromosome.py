#!/usr/bin/env python3
"""One chromosome-length artificial decode with the unchanged max automaton."""
from bisect import bisect_left, bisect_right
import json
import math
import os
from pathlib import Path
import resource
import time

import numpy as np
import torch
import manifest as M
import profile_emission_lattice_synthetic as S
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.full_grammar_viterbi import decode_full_grammar_viterbi, _forward, _trace

OUT = Path('outputs/M28-FULL-GRAMMAR-VITERBI-CHROMOSOME')
LENGTH = 30427671
CAP, CHECKPOINT = 128, 8192
OUTPUT_LIMIT = 64 * 1024 * 1024
CONFIG = dict(scope='single_artificial_contig_at_longest_original_TRAIN_chromosome_length',
    length=LENGTH, DNA_seed=0, potential_seed=41,
    DNA='IID_uniform_ACGT_not_research_FASTA',
    potentials='independent_FP32_integers_minus3_to3_per_original_raw_grid_window',
    width=24576, stride=12288, shifted_final_full_window=True, strand='+',
    ownership='nearest_containing_genomic_window_center_per_base_lower_start_ties',
    backend='unchanged_checkpointed_complete_boundary_state_automaton',
    max_exons=CAP, checkpoint_bp=CHECKPOINT, role_weights=1, null=0,
    score_accumulation='FP64_per_base_forward_order',
    ties='background_continuation_then_state_table_arrival',
    maximum_new_output_bytes=OUTPUT_LIMIT, optimizer_steps=0)


def save_json(name, value):
    encoded = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    if sum(p.stat().st_size for p in OUT.iterdir()) + len(encoded) > OUTPUT_LIMIT:
        raise ValueError('Declared artificial output limit exceeded')
    with (OUT/name).open('xb') as handle:
        handle.write(encoded)


def main():
    torch.set_num_threads(2)
    began = time.perf_counter(); job = os.environ['SLURM_JOB_ID']
    source_root = Path('outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-GLM-PREFLIGHT')
    source = json.loads((source_root/'summary.json').read_text())
    lengths = {species: row['TRAIN_sequence_lengths']
        for species,row in source['original_TRAIN_DNA_sources'].items()}
    if (not source['completed'] or source['job_id'] != '13429660'
            or (source_root/'JOBID').read_text().strip() != '13429660'
            or (source_root/'STATUS').read_text().strip() != 'COMPLETED'
            or set(lengths) != {'arabidopsis_thaliana','oryza_sativa'}
            or sum(len(v) for v in lengths.values()) != 6
            or lengths['arabidopsis_thaliana']['NC_003070.9'] != LENGTH
            or max(n for v in lengths.values() for n in v.values()) != LENGTH):
        raise ValueError('Completed original TRAIN length metadata changed')
    for case, gain in ((0,0.),(1,21526.)):
        root = Path('outputs/M28-FULL-GRAMMAR-VITERBI-SCALE')/f'case_{case}'
        previous = json.loads((root/'summary.json').read_text())
        if (not previous['completed'] or previous['job_id'] != f'13429908_{case}'
                or (root/'STATUS').read_text().strip() != 'COMPLETED'
                or (root/'JOBID').read_text().strip() != previous['job_id']
                or previous['run']['gain'] != gain
                or not previous['run']['optimal_integer_score_matches_independent_Python_reference']):
            raise ValueError('Actual completed fixed 49kb independent optimum required')
    if (M.WIDTH,M.STRIDE) != (24576,12288):
        raise ValueError('Original raw grid changed')
    started = dict(job_id=job, configuration=CONFIG, started_utc=S.utc_now(),
        length_source_job_id='13429660', length_source_seqid='NC_003070.9',
        length_metadata_only=True, previous_max_representation_jobs=['13429852','13429908_0','13429908_1'])
    save_json('started.json',started); print(json.dumps(started,allow_nan=False),flush=True)
    generator = torch.Generator().manual_seed(0)
    codes = torch.randint(0,4,(LENGTH,),generator=generator)
    sequence = np.frombuffer(b'ACGT',dtype=np.uint8)[codes.numpy()].tobytes().decode('ascii')
    del codes
    locations = [(a,min(a+M.WIDTH,LENGTH)) for a in M.starts(LENGTH)]
    potential_generator = torch.Generator().manual_seed(41)

    def records():
        for i,(a,b) in enumerate(locations):
            values = torch.randint(-3,4,(b-a,15),generator=potential_generator).float()
            yield dict(window_index=i,sequence=sequence[a:b],local_potentials=values)

    joined = stitch_contig_emissions(sequence,locations,'+',records())
    logits = joined['local_potentials']; upper = 0.
    for a in range(0,LENGTH,65536):
        chunk = logits[a:a+65536].double()
        upper += float((chunk-chunk[:,:1]).max(1).values.sum())
    del chunk
    preparation_seconds = time.perf_counter()-began
    print(json.dumps(dict(stage='chromosome_decode_started',utc=S.utc_now(),
        preparation_seconds=preparation_seconds,windows=len(locations),
        pointwise_upper_bound=upper,input_logit_bytes=logits.numel()*logits.element_size())),flush=True)
    began = time.perf_counter()
    result = decode_full_grammar_viterbi(sequence,logits,CAP,CHECKPOINT)
    decode_seconds = time.perf_counter()-began
    print(json.dumps(dict(stage='chromosome_decode_finished',utc=S.utc_now(),
        decode_seconds_including_cold_JIT=decode_seconds,selected_count=len(result['chains']),
        gain=result['gain'],counts=result['counts'])),flush=True)
    began = time.perf_counter(); previous_end=0; total=0.; maximum_exons=0; maximum_span=0
    crossing_checkpoint=0; crossing_context=0; written=0
    occupied=sum(p.stat().st_size for p in OUT.iterdir())
    cuts = [a for _,a,_ in joined['genomic_context_segments']][1:]
    if result['reference_used'] or not 0<=result['gain']<=upper:
        raise ValueError('Reference-free whole-contig score and pointwise bound required')
    with (OUT/'selected_chains.jsonl').open('xb') as handle:
        for chain,gain in zip(result['chains'],result['chain_gains']):
            a,b=chain[0][0],chain[-1][1]
            if a<previous_end: raise ValueError('Selected whole-contig genes overlap')
            local=tuple((x-a,y-a) for x,y in chain)
            # The existing independent verifier now sees only this selected span:
            # its total role allocation is bounded by disjoint gene spans, not N per gene.
            S.verify_selected(sequence[a:b],logits[a:b],
                dict(chains=(local,),chain_gains=(gain,),gain=gain,reference_used=False),CAP)
            encoded=(json.dumps(dict(chain=chain,gain=gain),allow_nan=False)+'\n').encode()
            written+=len(encoded)
            if occupied+written>OUTPUT_LIMIT:
                raise ValueError('Declared selected-chain output limit exceeded')
            handle.write(encoded); total+=gain; previous_end=b
            maximum_exons=max(maximum_exons,len(chain)); maximum_span=max(maximum_span,b-a)
            crossing_checkpoint+=a//CHECKPOINT!=(b-1)//CHECKPOINT
            crossing_context+=bisect_left(cuts,b)>bisect_right(cuts,a)
    if total!=result['gain'] or len(result['chains'])!=len(result['chain_gains']):
        raise ValueError('All complete selected-chain scores must equal the whole-path score')
    counts=result['counts']; states=counts['state_count']
    if (counts['checkpoint_bytes']!=(math.ceil(LENGTH/CHECKPOINT)+1)*states*8
            or counts['maximum_block_predecessor_bytes']!=CHECKPOINT*states*4
            or counts['selected_role_path_bytes']!=LENGTH
            or counts['input_logit_bytes']!=LENGTH*15*4
            or counts['forward_and_recomputed_bases']!=2*LENGTH
            or not _forward.nopython_signatures or not _trace.nopython_signatures):
        raise ValueError('Actual chromosome arrays and two kernels differ from fixed representation')
    run=dict(preparation_seconds=preparation_seconds,decode_seconds_including_cold_JIT=decode_seconds,
        selected_verification_seconds=time.perf_counter()-began,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        gain=result['gain'],pointwise_upper_bound=upper,selected_count=len(result['chains']),
        maximum_selected_exons=maximum_exons,maximum_selected_gene_span=maximum_span,
        selected_genes_crossing_checkpoint_boundary=crossing_checkpoint,
        selected_genes_crossing_owned_context_boundary=crossing_context,counts=counts,
        windows=len(locations),first_window=locations[0],last_window=locations[-1],
        motif_counts_full_substrings=dict(start=sequence.count('ATG'),donor=sequence.count('GT')+sequence.count('GC'),
            acceptor=sequence.count('AG'),stop=sum(sequence.count(v) for v in S.STOP)),
        forward_nopython_signatures=len(_forward.nopython_signatures),trace_nopython_signatures=len(_trace.nopython_signatures),
        all_selected_ORF_carry_contiguous_terminal_splice_cap_nonoverlap_and_integer_scores_recomputed=True,
        chromosome_length_independent_exhaustive_optimum_available=False)
    summary=dict(**started,completed=True,completed_utc=S.utc_now(),run=run,
        research_DNA_or_labels_opened=False,original_TRAIN_length_metadata_read=True,
        model_or_backbone_forward=False,checkpoint_read=False,optimizer_steps=0,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,
        stopped_A43_endpoint_prototype_used=False,A44_execution_chain_reopened=False,
        interpretation='chromosome_length_artificial_exact_representation_resource_not_research_annotator_performance')
    save_json('summary.json',summary); print(json.dumps(summary,allow_nan=False),flush=True)


if __name__=='__main__':
    main()
