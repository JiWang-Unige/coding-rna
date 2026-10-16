#!/usr/bin/env python3
"""A21 raw-FASTA caller; actual TRAIN PASS and a separate DEV protocol required."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import time
import torch
import run_junction_chain_set as J
from src.m28.training import one_hot
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

C = J.R.C
ROOT = J.ROOT
M = C.I.M
DEV_PROTOCOL = 'reports/M28-METHOD-RESTART/a21_junction_chain_dev_protocol.md'
SCORE_KIND = 'A21_junction_chain_logit_minus_window_null'


def validate_train_preparation(fitted, prepared, checked):
    J.validate_fit(fitted)
    if (fitted['quality_trainable_parameters'] != 248068 or
        fitted['generator_trainable_parameters'] != 0 or
        [e['steps'] for e in fitted['epochs']] != [1536]*3):
        raise ValueError('Fixed A21 architecture and original 4608-step exposure required')
    for receipt in (prepared, checked):
        if (receipt['completed'] is not True or receipt['configuration'] != J.CONFIG or
            receipt['protocol'] != J.PROTOCOL or receipt['windows'] != 1491 or
            receipt['candidate_window_pairs'] != 16321 or receipt['optimizer_steps_this_run'] != 0 or
            receipt['DEV_forward_used'] is not False or receipt['test_setaria_used'] is not False):
            raise ValueError('Complete original natural TRAIN domain required')
    selected = prepared['natural_selection_per_species']
    if set(selected) != set(C.R.ALLOW) or checked['natural_selection_per_species'] != selected:
        raise ValueError('Both species and reconciled TRAIN selections required')
    checks = J.gate(selected)
    if (prepared['checkpoint'] != fitted['final_checkpoint'] or
        prepared['predeclared_TRAIN_preparation_gate'] != checks or
        checked['original_frozen_TRAIN_gate'] != checks or
        prepared['TRAIN_preparation_supported'] is not True or
        checked['TRAIN_preparation_supported'] is not True or not all(checks.values())):
        raise ValueError('The original five-condition A21 TRAIN preparation gate must pass')
    if (prepared['threshold'] != 0. or prepared['reference_injected_candidates'] != 0 or
        any(prepared[k] is not False for k in ('new_candidate_search', 'new_backbone_forward', 'threshold_scan')) or
        any(prepared[k] is not True for k in ('natural_candidate_geometry_labels_and_origins_unchanged',
                                            'same_fit_immutable_hidden_used')) or
        any(checked[k] is not True for k in ('actual_4608_steps_draw_weights_and_original_natural_domain_verified',
            'actual_last_head_shapes_finite_tensors_components_selection_profiles_and_gate_verified',
            'cache_zero_update_equivalence_verified', 'original_unknown_semantics_unchanged')) or
        checked['new_head_or_backbone_inference'] is not False):
        raise ValueError('Unchanged candidates, operating point and actual result reconciliation required')


@torch.inference_mode()
def quality_record(generator, head, hidden, sequence, index, ident):
    chains, sources = C.free_pool(generator, hidden, sequence)
    result = head(hidden, chains, len(sequence))
    components = {k: result[k].float().cpu().tolist() for k in
                  ('exon_additive', 'junction_additive', 'nonadditive')}
    null = float(result['null']); scores = result['logits'].float().cpu().tolist()
    if any(not math.isfinite(v) for v in [null]+scores+[v for values in components.values() for v in values]):
        raise ValueError('Nonfinite A21 chain score component')
    return dict(window_index=index, window_id=ident, chains=chains, score=scores,
        decoder_sources=sources, score_components=components, window_null=null,
        score_kind=SCORE_KIND, quality_logit_threshold=0., reference_used=False,
        counts=dict(free_candidates=len(chains), nonpositive_candidates=sum(v <= 0 for v in scores)))


@torch.inference_mode()
def annotate_contig(generator, head, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    if not sequence:
        raise ValueError('Empty contig')
    length = len(sequence)
    locations = [[a, min(a+M.WIDTH, length)] for a in M.starts(length)]
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
                free = quality_record(generator, head, hidden, oriented, index, ident)
                handle.write(json.dumps(free, allow_nan=False)+'\n'); counts.update(free['counts'])
                if (index+1) % 32 == 0:
                    handle.flush(); print(json.dumps(dict(species=species, strand=strand,
                        windows=index+1, total=len(locations))), flush=True)
                del features, dna, hidden, free
        assembled = assemble_predictions(C.I.saved_records(path), locations, strand)
        for index, row in enumerate(assembled['predictions']):
            ident = 'M28_A21_'+str(ordinal)+'_'+strand+'_'+str(index)
            transcripts.append(as_transcript(row['chain'], seqid, strand, ident, row['score']))
        scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
            windows=locations, file=path.name, quality_counts=dict(counts),
            assembly_counts=assembled['counts'], seconds=time.perf_counter()-started))
    return transcripts, scopes


def main(out, fasta=None):
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (J.FIT, J.READOUT, J.CHECK)):
        raise ValueError('Terminal A21 fit, TRAIN readout and CPU reconciliation required')
    fitted, prepared, checked = [json.loads((p/'summary.json').read_text()) for p in (J.FIT, J.READOUT, J.CHECK)]
    validate_train_preparation(fitted, prepared, checked)
    if ((J.FIT/'JOBID').read_text().strip() != (J.READOUT/'JOBID').read_text().strip() or
        checked['source_fit_readout_job'] != (J.FIT/'JOBID').read_text().strip()):
        raise ValueError('CPU reconciliation must describe the actual fit/readout job')
    # This is a scientific contract prerequisite, not evidence of human approval.
    if not (ROOT/DEV_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed A21 raw-FASTA DEV protocol before native execution')
    if out.exists():
        raise FileExistsError('A21 raw inference requires a new output directory')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    J.validate_fit(state, state=True)
    head = J.new_head('cuda').eval().requires_grad_(False)
    head.load_state_dict(state['head'], strict=True); del state
    generator = C.load_generator(); backbone, tokenizer = C.I.load_native()
    out.mkdir(parents=True, exist_ok=False)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (species, seqid, sequence) in enumerate(C.I.input_contigs(fasta)):
        txs, ss = annotate_contig(generator, head, backbone, tokenizer, species, seqid,
                                  sequence, ordinal, out, 'cuda')
        transcripts.extend(txs); scopes.extend(ss); lengths[seqid] = len(sequence)
    windows = sum(len(s['windows']) for s in scopes)
    if fasta is None and windows != 8690:
        raise ValueError('Incomplete natural DEV grid')
    write_gff3(out/'A21.gff3', transcripts, lengths, 'M28_A21')
    C.R.write_summary(out, dict(arm='A21', completed=True, inference_complete=True,
        configuration=J.CONFIG, protocol=DEV_PROTOCOL, TRAIN_protocol=J.PROTOCOL,
        TRAIN_preparation_gate_passed=True, checkpoint=fitted['final_checkpoint'], optimizer_steps=4608,
        generator_checkpoint=C.GENERATOR, feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s: str(ROOT/'data/m1_screen'/s/'genome.fa') for s in C.I.DEV_LENGTHS},
        reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
        DEV_used_for_selection=False, test_setaria_used=False, windows=windows, scopes=scopes, lengths=lengths,
        primary_CDS_chain_predictions=len(transcripts), wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='Raw FASTA annotation is not an accuracy or generalization result before separate evaluation.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev', action='store_true'); group.add_argument('--fasta', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(), args.fasta.resolve() if args.fasta else None)
