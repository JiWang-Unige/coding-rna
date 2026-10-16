#!/usr/bin/env python3
"""A36 raw caller: unchanged A3 natural proposals, separately fitted emissions."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import time
import torch
import readout_structured_emission_train as D
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3
from src.m28.training import one_hot

ROOT, C, M = D.ROOT, D.R.C, D.R.C.I.M
READOUT_JOB, CHECK_JOB = '13408449', '13408451'
DEV_PROTOCOL = 'reports/M28-METHOD-RESTART/a36_structured_emission_dev_protocol.md'
SCORE_KIND = 'A36_fixed_unit_structured_role_emission_gain'
CPU_PROOFS = (
    'actual_original_training_draw_weights_and_4608_updates_verified',
    'actual_three_checkpoint_tensors_and_optimizer_exposure_verified',
    'actual_saved_natural_identity_and_complete_logits_verified',
    'actual_saved_full_window_log_probability_scores_verified',
    'actual_saved_labels_selections_stronger_baseline_and_five_conditions_recomputed',
)


def validate_preparation(fitted, prepared, checked):
    D.validate_fit(fitted)
    for report in (prepared, checked):
        if (report['completed'] is not True or report['configuration'] != D.CONFIG or
                report['protocol'] != D.PROTOCOL or report['windows'] != 1491 or
                report['natural_candidates'] != 16321 or report['reference_injected_predictions'] != 0 or
                report['optimizer_steps_this_run'] != 0 or any(report[k] is not False for k in
                ('new_candidate_search', 'new_R4_labelling', 'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
            raise ValueError('Complete unchanged original TRAIN-only A36 preparation required')
    if (prepared['source_fit_job'] != D.FIT_JOB or prepared['checkpoint'] != fitted['final_checkpoint'] or
            checked['source_jobs'] != dict(fit=D.FIT_JOB, readout=READOUT_JOB) or
            any(checked[k] is not True for k in CPU_PROOFS) or checked['model_or_backbone_forward'] is not False or
            prepared['label_catalog_opened'] is not False or prepared['new_backbone_forward'] is not False or
            any(prepared[k] is not True for k in ('natural_score_file_closed_before_saved_label_arrays_open',
            'original_natural_geometry_order_sources_unchanged', 'original_saved_labels_reused',
            'model_weights_unchanged', 'all_parameters_frozen_and_gradients_absent'))):
        raise ValueError('Actual matching saved score, model and CPU reconciliation required')
    counts, refs = prepared['selected_counts_per_species'], prepared['eligible_reference_pairs']
    if (set(counts) != set(C.R.ALLOW) or set(refs) != set(C.R.ALLOW) or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs or
            sum(v['eligible'] for v in refs.values()) != 2916):
        raise ValueError('Both species, original denominators and paired stronger baseline required')
    gates = D.gate(counts)
    if (prepared['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            prepared['TRAIN_preparation_supported'] is not True or
            checked['TRAIN_preparation_supported'] is not True or not all(gates.values())):
        raise ValueError('Original five-condition A36 TRAIN gate must actually pass')


def readiness():
    D.F.A.require_jobs([D.FIT_JOB, READOUT_JOB, CHECK_JOB])
    fitted = D.R.completed(D.F.OUT, D.FIT_JOB)
    prepared = D.R.completed(D.OUT, READOUT_JOB); checked = D.R.completed(D.CHECK, CHECK_JOB)
    validate_preparation(fitted, prepared, checked)
    return fitted


def load_models(fitted):
    generator = C.load_generator()
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    D.validate_fit(state, state=True)
    emission = D.F.A.new_model('cuda').eval().requires_grad_(False)
    emission.load_state_dict(state['model'], strict=True)
    return generator, emission


@torch.inference_mode()
def local_record(generator, hidden, logits, sequence, index, ident):
    chains, sources = C.free_pool(generator, hidden, sequence)
    scores = D.local_path_scores(logits, chains, len(sequence)).cpu().tolist()
    return dict(window_index=index, window_id=ident, chains=chains, score=scores,
        decoder_sources=sources, score_kind=SCORE_KIND, variant='A36', score_threshold=0.,
        reference_used=False, counts=dict(free_candidates=len(chains),
        nonpositive_candidates=sum(v <= 0 for v in scores)))


@torch.inference_mode()
def annotate_contig(generator, emission, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
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
                logits = emission(features, dna, [b-a])['segmentation_logits'][0]
                if logits.dtype != torch.float32:
                    raise ValueError('Original FP32 emission precision required')
                ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                record = local_record(generator, hidden, logits, oriented, index, ident)
                handle.write(json.dumps(record, allow_nan=False)+'\n'); counts.update(record['counts'])
                if (index+1) % 32 == 0:
                    handle.flush(); print(json.dumps(dict(variant='A36', species=species,
                        strand=strand, windows=index+1, total=len(locations))), flush=True)
                del features, dna, hidden, logits, record
        assembled = assemble_predictions(C.I.saved_records(path), locations, strand)
        for index, row in enumerate(assembled['predictions']):
            ident = 'M28_A36_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
            transcripts.append(as_transcript(row['chain'], seqid, strand, ident, row['score']))
        scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
            windows=locations, file=path.name, candidate_counts=dict(counts),
            assembly_counts=assembled['counts'], seconds=time.perf_counter()-started))
    return transcripts, scopes


@torch.inference_mode()
def main(out, fasta=None):
    fitted = readiness()
    # A scientific contract, never evidence of human permission to access new data.
    if not (ROOT/DEV_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed A36 raw-FASTA DEV protocol before native execution')
    if out.exists():
        raise FileExistsError('A36 raw inference requires a new output directory')
    generator, emission = load_models(fitted); backbone, tokenizer = C.I.load_native()
    models = dict(generator=generator, emission=emission)
    before = {a: {k: v.detach().clone() for k, v in m.state_dict().items()} for a, m in models.items()}
    out.mkdir(parents=True, exist_ok=False)
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (species, seqid, sequence) in enumerate(C.I.input_contigs(fasta)):
        if seqid in lengths:
            raise ValueError('Raw FASTA contig identifiers must be unique for GFF coordinates')
        txs, local = annotate_contig(generator, emission, backbone, tokenizer, species, seqid,
                                     sequence, ordinal, out, 'cuda')
        transcripts.extend(txs); scopes.extend(local); lengths[seqid] = len(sequence)
    windows = sum(len(s['windows']) for s in scopes)
    if not lengths:
        raise ValueError('Raw FASTA must contain at least one nonempty contig')
    if fasta is None and windows != 8690:
        raise ValueError('Incomplete natural DEV grid')
    for arm, model in models.items():
        if (any(not torch.equal(v, before[arm][k]) for k, v in model.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen raw caller altered parameters or received gradients')
    write_gff3(out/'A36.gff3', transcripts, lengths, 'M28_A36')
    C.R.write_summary(out, dict(arm='A36', completed=True, inference_complete=True,
        configuration=D.CONFIG, protocol=DEV_PROTOCOL, TRAIN_protocol=D.PROTOCOL,
        fit_protocol=D.F.PROTOCOL, TRAIN_preparation_gate_passed=True,
        source_TRAIN_jobs=dict(fit=D.FIT_JOB, readout=READOUT_JOB, CPU=CHECK_JOB),
        generator_checkpoint=C.GENERATOR, emission_checkpoint=fitted['final_checkpoint'],
        generator_optimizer_steps=4608, emission_optimizer_steps=4608,
        model_parameters={a: sum(p.numel() for p in m.parameters()) for a, m in models.items()},
        feature_revision=C.I.F.REVISION, feature_spec=C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s: str(ROOT/'data/m1_screen'/s/'genome.fa') for s in C.I.DEV_LENGTHS},
        score_kind=SCORE_KIND, reference_used=False, labels_or_feature_cache_used=False,
        GEN_encoded_once_per_oriented_window=True, natural_proposal_generator='unchanged_A3_greedy_plus_beam4',
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0, DEV_used_for_selection=False, test_setaria_used=False,
        A33_outputs_opened=False, windows=windows, scopes=scopes, lengths=lengths,
        primary_CDS_chain_predictions=len(transcripts), wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Raw primary CDS caller with additional learned representation and training exposure; '
            'not matched-capacity ablation, accuracy, causal novelty or independent generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev', action='store_true'); group.add_argument('--fasta', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(), args.fasta.resolve() if args.fasta else None)
