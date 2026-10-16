#!/usr/bin/env python3
"""Real frozen encoder/caller on artificial short contigs; no biological evaluation."""
import json
from pathlib import Path
import time
import torch
import infer_local_path_caller as I
import evaluate_chain_quality as E
from src.m28.export import write_gff3

OUT = I.ROOT/'outputs/M28-TRANSFER-SHORT-RAW-PREFLIGHT'


@torch.inference_mode()
def main():
    prepared, checked = I.R.completed(I.R.READOUT), I.R.completed(I.R.CHECK)
    for arm in ('A3', 'A23'):
        I.validate_train_preparation(prepared, checked, arm)
    if checked['source_jobs'] != {'readout': (I.R.READOUT/'JOBID').read_text().strip()}:
        raise ValueError('Actual original A31 TRAIN readout identity required')
    fixtures = {
        'artificial_332': 'ATG'+'CAA'*108+'TAA'+'cc',
        'artificial_7536': 'C'*6140+'ATG'+'CAA'*463+'TAA'+'C',
    }
    lengths = {name: len(sequence) for name, sequence in fixtures.items()}
    if list(lengths.values()) != [332, 7536]:
        raise ValueError('Fixed short-contig engineering lengths changed')
    fasta = OUT/'artificial.fa'
    with fasta.open('x') as handle:
        for name, sequence in fixtures.items():
            handle.write('>'+name+' artificial DNA, not an organism\n'+sequence+'\n')
    contigs = list(I.C.I.input_contigs(fasta))
    if [(name, sequence) for _, name, sequence in contigs] != [
            (name, sequence.upper()) for name, sequence in fixtures.items()]:
        raise ValueError('Explicit multi-record FASTA identity/normalization changed')
    began = time.perf_counter()
    backbone, tokenizer = I.C.I.load_native()
    torch.cuda.reset_peak_memory_stats()
    variants = {}
    for arm in ('A3', 'A23'):
        model = I.load_generator(arm)
        before = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        observed = []

        def observe(module, args, output):
            features, dna, valid = args
            if (tuple(features.shape) != (1, 4096, 2048) or
                    tuple(dna.shape) != (1, I.M.WIDTH, 5) or len(valid) != 1 or
                    not bool(torch.isfinite(output['features']).all())):
                raise ValueError('Real padded encoder/head tensor interface failed')
            observed.append(dict(valid_bases=valid[0], feature_shape=list(features.shape),
                dna_shape=list(dna.shape), hidden_shape=list(output['features'].shape),
                feature_dtype=str(features.dtype), hidden_finite=True))

        hook = model.register_forward_hook(observe)
        directory = OUT/arm
        directory.mkdir()
        transcripts, scopes = [], []
        try:
            for ordinal, (species, seqid, sequence) in enumerate(contigs):
                txs, local = I.annotate_contig(model, backbone, tokenizer, species, seqid,
                    sequence, ordinal, directory, 'cuda', arm)
                transcripts.extend(txs)
                scopes.extend(local)
        finally:
            hook.remove()
        if [r['valid_bases'] for r in observed] != [332, 332, 7536, 7536]:
            raise ValueError('Both independent contigs and strands must reach the real head')
        rows = 0
        for scope in scopes:
            length = lengths[scope['seqid']]
            if scope['windows'] != [[0, length]]:
                raise ValueError('Short contig omitted or padded bases exposed as input coordinates')
            records = list(I.C.I.saved_records(directory/scope['file']))
            if (len(records) != 1 or records[0]['window_index'] != 0 or
                    records[0]['window_id'] != 'user_input|'+scope['seqid']+'|0|'+scope['strand']):
                raise ValueError('Saved short-contig identity or scope mismatch')
            for chain in records[0]['chains']:
                if any(not 0 <= a < b <= length for a, b in chain):
                    raise ValueError('Natural chain extends beyond the real contig')
            rows += len(records)
        if sum(s['assembly_counts']['selected_unique'] for s in scopes) != len(transcripts):
            raise ValueError('Whole-contig assembly and GFF population mismatch')
        if (any(p.requires_grad or p.grad is not None for p in model.parameters()) or
                any(not torch.equal(before[k], v.detach().cpu()) for k, v in model.state_dict().items())):
            raise ValueError('Frozen generator changed during engineering inference')
        gff = directory/('A31_'+arm+'.gff3')
        write_gff3(gff, transcripts, lengths, 'M28_A31_'+arm)
        parsed = E.A.E.S.E.primary_transcripts(E.A.E.S.E.parse_annotation(gff, lengths), complete_only=False)
        E.A.validate_gff_roundtrip(parsed, transcripts)
        variants[arm] = dict(scopes=scopes, saved_rows=rows, real_head_inputs=observed,
            synthetic_prediction_count=len(transcripts), gff=str(gff.relative_to(I.ROOT)),
            original_generator_state_unchanged=True, gff_ids_coordinates_phase_roundtrip=True)
        del model, before
    report = dict(experiment='M28-TRANSFER-SHORT-RAW-PREFLIGHT', completed=True,
        input_kind='artificial_DNA_only', input_lengths=lengths, variants=variants,
        encoder_model=I.C.I.F.MODEL, encoder_revision=I.C.I.F.REVISION,
        feature_spec=I.C.R.FEATURE_SPEC, caller_configuration=I.R.CONFIG,
        original_A31_TRAIN_preparation_identity_verified=True,
        real_native_encoder=True, new_training_updates=0, new_parameters=0,
        biological_genome_or_annotation_opened=False, feature_cache_used=False,
        DEV_outputs_opened=False, test_setaria_used=False, independent_panel_selected=False,
        existing_inference_code_or_protocol_changed=False,
        wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Engineering readiness at two artificial lengths, not accuracy, throughput, generalization or model selection.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
