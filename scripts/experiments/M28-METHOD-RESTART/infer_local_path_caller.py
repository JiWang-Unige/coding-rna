#!/usr/bin/env python3
"""A31 raw-FASTA interface; actual variant TRAIN PASS and separate DEV contract required."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import time
import torch
import run_local_path_caller as R
from src.m28.local_path_score import local_path_scores
from src.m28.training import one_hot
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

C, ROOT = R.C, R.ROOT
M = C.I.M
DEV_PROTOCOL = 'reports/M28-METHOD-RESTART/a31_local_path_caller_dev_protocol.md'
SCORE_KIND = 'A31_fixed_local_role_emission_gain'


def validate_train_preparation(prepared, checked, generator):
    arm = 'A31_'+generator
    for receipt in (prepared, checked):
        if (receipt['completed'] is not True or receipt['arm'] != 'A31' or
                receipt['configuration'] != R.CONFIG or receipt['protocol'] != R.PROTOCOL or
                receipt['windows'] != 1491 or receipt['optimizer_steps_this_run'] != 0 or
                receipt['reference_injected_candidates'] != 0 or
                any(receipt[k] is not False for k in
                    ('new_candidate_search', 'new_backbone_forward', 'DEV_forward_used', 'test_setaria_used'))):
            raise ValueError('Matching frozen TRAIN-only A31 readout and CPU reconciliation required')
    counts, references = prepared['selected_counts_per_species'], prepared['eligible_reference_pairs']
    if (set(counts) != set(C.R.ALLOW) or set(references) != set(C.R.ALLOW) or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != references):
        raise ValueError('Both species and reconciled reference denominators required')
    gates = R.gate(counts, references)
    if (prepared['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            prepared['TRAIN_preparation_supported'][arm] is not True or
            checked['TRAIN_preparation_supported'][arm] is not True or not all(gates[arm].values())):
        raise ValueError('Selected variant must pass the original five-condition A31 TRAIN gate')
    if (any(prepared[k] is not True for k in
            ('natural_score_file_closed_before_catalog_open', 'original_greedy_source_replay_verified',
             'model_weights_unchanged', 'all_parameters_frozen_and_gradients_absent', 'original_natural_pools_unchanged')) or
            any(checked[k] is not True for k in
            ('actual_saved_natural_pool_and_greedy_subset_identity_verified',
             'actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed',
             'actual_five_fixed_gate_integer_arithmetic_verified',
             'original_last_checkpoint_tensors_finite_and_parameter_counts_verified')) or
            checked['new_head_or_backbone_forward'] is not False):
        raise ValueError('Actual frozen model, natural pool and metric reconciliation required')


def load_generator(generator):
    if generator == 'A3':
        return C.load_generator()
    fitted = R.completed(R.R23.FIT, 13391608); R.R23.validate_fit(fitted)
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R.R23.validate_fit(state, state=True)
    model = R.R23.new_model().eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True)
    return model


@torch.inference_mode()
def local_record(generator, hidden, sequence, index, ident, arm):
    chains, sources = C.free_pool(generator, hidden, sequence)
    scores = local_path_scores(generator.segmentation(hidden), chains, len(sequence)).cpu().tolist()
    if not all(math.isfinite(v) for v in scores):
        raise ValueError('Nonfinite fixed local emission score')
    # Keep the entire natural pool, including zero and negative gains.
    return dict(window_index=index, window_id=ident, chains=chains, score=scores,
        decoder_sources=sources, score_kind=SCORE_KIND, variant='A31_'+arm,
        score_threshold=0., reference_used=False,
        counts=dict(free_candidates=len(chains), nonpositive_candidates=sum(v <= 0 for v in scores)))


@torch.inference_mode()
def annotate_contig(generator, backbone, tokenizer, species, seqid, sequence, ordinal, out, device, arm):
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    if not sequence:
        raise ValueError('Empty contig')
    locations = [[a, min(a+M.WIDTH, len(sequence))] for a in M.starts(len(sequence))]
    transcripts, scopes = [], []
    for strand in ('+', '-'):
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand == '+' else 'm')+'.jsonl')
        started = time.perf_counter(); counts = Counter()
        with path.open('x') as handle:
            for index, (a, b) in enumerate(locations):
                raw = sequence[a:b]; oriented = raw if strand == '+' else M.rc(raw)
                features = C.I.encode_window(backbone, tokenizer, oriented, device)
                dna = one_hot(oriented, M.WIDTH)[None].to(device)
                hidden = generator(features, dna, [b-a])['features'][0]
                ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                record = local_record(generator, hidden, oriented, index, ident, arm)
                handle.write(json.dumps(record, allow_nan=False)+'\n'); counts.update(record['counts'])
                if (index+1) % 32 == 0:
                    handle.flush(); print(json.dumps(dict(variant='A31_'+arm, species=species,
                        strand=strand, windows=index+1, total=len(locations))), flush=True)
                del features, dna, hidden, record
        assembled = assemble_predictions(C.I.saved_records(path), locations, strand)
        for index, row in enumerate(assembled['predictions']):
            ident = 'M28_A31_'+arm+'_'+str(ordinal)+'_'+strand+'_'+str(index)
            transcripts.append(as_transcript(row['chain'], seqid, strand, ident, row['score']))
        scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
            windows=locations, file=path.name, candidate_counts=dict(counts),
            assembly_counts=assembled['counts'], seconds=time.perf_counter()-started))
    return transcripts, scopes


def main(out, generator, fasta=None):
    prepared = R.completed(R.READOUT); checked = R.completed(R.CHECK)
    validate_train_preparation(prepared, checked, generator)
    if checked['source_jobs'] != {'readout': (R.READOUT/'JOBID').read_text().strip()}:
        raise ValueError('CPU reconciliation must describe the actual readout job')
    # A scientific contract prerequisite, not evidence of human approval.
    if not (ROOT/DEV_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed A31 raw-FASTA DEV protocol before native execution')
    if out.exists():
        raise FileExistsError('A31 raw inference requires a new output directory')
    model = load_generator(generator); backbone, tokenizer = C.I.load_native()
    out.mkdir(parents=True, exist_ok=False)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (species, seqid, sequence) in enumerate(C.I.input_contigs(fasta)):
        txs, local = annotate_contig(model, backbone, tokenizer, species, seqid, sequence,
                                     ordinal, out, 'cuda', generator)
        transcripts.extend(txs); scopes.extend(local); lengths[seqid] = len(sequence)
    windows = sum(len(s['windows']) for s in scopes)
    if fasta is None and windows != 8690:
        raise ValueError('Incomplete natural DEV grid')
    write_gff3(out/('A31_'+generator+'.gff3'), transcripts, lengths, 'M28_A31_'+generator)
    C.R.write_summary(out, dict(arm='A31', variant='A31_'+generator, completed=True,
        inference_complete=True, configuration=R.CONFIG, protocol=DEV_PROTOCOL, TRAIN_protocol=R.PROTOCOL,
        TRAIN_preparation_gate_passed=True, generator_checkpoint=R.CONFIG['generators'][generator],
        generator_optimizer_steps=4608, feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s: str(ROOT/'data/m1_screen'/s/'genome.fa') for s in C.I.DEV_LENGTHS},
        score_kind=SCORE_KIND, reference_used=False, labels_or_feature_cache_used=False,
        optimizer_steps_this_run=0, new_parameters=0, DEV_used_for_selection=False,
        test_setaria_used=False, windows=windows, scopes=scopes, lengths=lengths,
        primary_CDS_chain_predictions=len(transcripts), wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Raw-FASTA annotation is not an accuracy or independent generalization result before evaluation.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev', action='store_true'); group.add_argument('--fasta', type=Path)
    parser.add_argument('--generator', choices=('A3', 'A23'), required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(), args.generator, args.fasta.resolve() if args.fasta else None)
