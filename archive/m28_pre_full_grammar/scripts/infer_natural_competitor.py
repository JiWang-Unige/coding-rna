#!/usr/bin/env python3
"""A43 raw interface preparation: one emission model, joint contig grammar."""
import argparse
import json
from pathlib import Path
import re
import time
import torch

import readout_natural_competitor_train as D
from src.m28.contig_emission_lattice import decode_contig_emission_lattice
from src.m28.export import as_transcript, write_gff3
from src.m28.training import one_hot

ROOT, NATIVE, M = D.ROOT, D.R.C.I, D.R.C.I.M
POTENTIAL_JOB, DECODE_JOB, READOUT_JOB, CHECK_JOB = '13424395', '13424400', '13424410', '13424413'
RAW_PROTOCOL = 'reports/M28-METHOD-RESTART/a43_natural_competitor_raw_protocol.md'
CONTEXT_POLICY = 'nearest_containing_genomic_window_center_per_base_lower_start_ties'
SCORE_KIND = 'A43_unit_role_gain_under_joint_base_owned_contig_context'
CPU_PROOFS = (
    'actual_original_draw_weights_union_NLL_and_4608_updates_verified',
    'actual_three_finite_A43_checkpoints_and_AdamW_exposure_verified',
    'actual_frozen_last_model_potential_identity_and_complete_32_shards_verified',
    'actual_new_DNA_ORF_motif_cap_span_and_full_window_FP64_log_probability_scores_verified',
    'actual_full_grammar_objective_not_below_same_model_original_natural_pool',
    'actual_frozen_R4_unknowns_stronger_baseline_and_complete_denominators_recomputed',
    'actual_per_species_counts_and_original_five_conditions_recomputed',
)


def validate_preparation(fitted, prepared, checked):
    D.validate_fit(fitted)
    source = dict(fit=D.FIT_JOB, potentials=POTENTIAL_JOB, decode=DECODE_JOB)
    for report in (prepared, checked):
        if (report['completed'] is not True or report['configuration'] != D.CONFIG or
                report['protocol'] != D.PROTOCOL or report['windows'] != 1491 or
                report['original_draws'] != 1536 or report['original_candidates'] != 16321 or
                report['reference_injected_predictions'] != 0 or report['optimizer_steps_this_run'] != 0 or
                report['new_geometry_labels_computed'] is not True or
                any(report[k] is not False for k in ('model_or_backbone_forward', 'new_candidate_search',
                    'R4_policy_or_catalog_changed', 'DEV_used', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened'))):
            raise ValueError('Actual matching original TRAIN-only A43 readout and CPU reconciliation required')
    if (prepared['source_jobs'] != source or checked['source_jobs'] != dict(source, readout=READOUT_JOB) or
            any(checked[k] is not True for k in CPU_PROOFS) or
            prepared['all_natural_search_and_score_files_closed_before_labels'] is not True or
            prepared['original_R4_labels_exact'] is not True or checked['candidate_window_pairs'] != 19189):
        raise ValueError('Actual A43 full-search sources and seven CPU proofs required')
    counts, refs = prepared['selected_counts_per_species'], prepared['eligible_reference_pairs']
    expected = {
        'arabidopsis_thaliana': (1531, 347, 315, 548, 2172, 854),
        'oryza_sativa': (478, 251, 132, 152, 744, 232),
    }
    if (set(counts) != set(expected) or set(refs) != set(expected) or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs):
        raise ValueError('Original stronger baseline and full species denominators required')
    for species, (correct, nonmatch, unknown, long, eligible, five) in expected.items():
        baseline = dict(selected=correct+nonmatch+unknown, selected_correct=correct,
            selected_known_nonmatch=nonmatch, selected_unknown=unknown, selected_five_or_more_correct=long)
        if counts[species]['A31_A3'] != baseline or refs[species] != dict(eligible=eligible, five_or_more=five):
            raise ValueError('Original stronger baseline and full species denominators required')
    gates = D.gate(counts)
    if (prepared['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            prepared['TRAIN_preparation_supported'] is not True or checked['TRAIN_preparation_supported'] is not True or
            not all(gates.values())):
        raise ValueError('Original five-condition A43 TRAIN gate must actually pass')


def readiness():
    require = D.F.A.B.D.B.F.A.require_jobs
    require([D.FIT_JOB, POTENTIAL_JOB, READOUT_JOB, CHECK_JOB] + [DECODE_JOB+'_'+str(i) for i in range(32)])
    fitted = D.R.completed(D.F.OUT, D.FIT_JOB)
    prepared = D.R.completed(D.OUT, READOUT_JOB); checked = D.R.completed(D.CHECK, CHECK_JOB)
    validate_preparation(fitted, prepared, checked)
    return fitted


def load_model(fitted, device):
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    D.validate_fit(state, state=True)
    model = D.F.A.A36.new_model(device).eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True)
    if sum(p.numel() for p in model.parameters()) != 366191:
        raise ValueError('Only the fixed A43 emission representation is allowed')
    return model


@torch.inference_mode()
def emission_records(model, backbone, tokenizer, sequence, locations, strand, device):
    """Only DNA and native encoding; no query model, pool, annotation or cache."""
    for index, (a, b) in enumerate(locations):
        raw = sequence[a:b]; oriented = raw if strand == '+' else M.rc(raw)
        features = NATIVE.encode_window(backbone, tokenizer, oriented, device)
        dna = one_hot(oriented, M.WIDTH)[None].to(device)
        logits = model(features, dna, [b-a])['segmentation_logits'][0][:b-a].cpu()
        yield dict(window_index=index, sequence=oriented, local_potentials=logits)
        del features, dna, logits


@torch.inference_mode()
def annotate_contig(model, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    if not sequence:
        raise ValueError('Empty contig')
    locations = [[a, min(a+M.WIDTH, len(sequence))] for a in M.starts(len(sequence))]
    transcripts, scopes = [], []
    for strand in ('+', '-'):
        began = time.perf_counter()
        records = emission_records(model, backbone, tokenizer, sequence, locations, strand, device)
        natural = decode_contig_emission_lattice(sequence, locations, strand, records, D.CONFIG['max_exons'])
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand == '+' else 'm')+'.json')
        with path.open('x') as handle:
            json.dump(dict(variant='A43', score_kind=SCORE_KIND, species=species, seqid=seqid, strand=strand,
                length=len(sequence), windows=locations, reference_used_in_scoring=False, **natural),
                handle, allow_nan=False); handle.write('\n')
        for index, (chain, gain) in enumerate(zip(natural['genomic_chains'], natural['chain_gains'])):
            ident = 'M28_A43_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
            transcripts.append(as_transcript(chain, seqid, strand, ident, gain))
        scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
            windows=locations, file=path.name, context_policy=natural['context_policy'],
            primary_CDS_chain_predictions=len(natural['chains']), seconds=time.perf_counter()-began))
    return transcripts, scopes


@torch.inference_mode()
def main(out, fasta):
    fitted = readiness()
    # A separate scientific/resource contract is required, not proof of human approval.
    if not (ROOT/RAW_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed A43 raw scientific/resource protocol before native execution')
    if fasta is None:
        raise ValueError('Explicit authorized FASTA required; no automatic DEV access')
    if out.exists():
        raise FileExistsError('A43 raw inference requires a new output directory')
    model = load_model(fitted, 'cuda'); backbone, tokenizer = NATIVE.load_native()
    backbone.eval().requires_grad_(False)
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    out.mkdir(parents=True, exist_ok=False); began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (species, seqid, sequence) in enumerate(NATIVE.input_contigs(fasta)):
        if seqid in lengths:
            raise ValueError('Raw FASTA contig identifiers must be unique for GFF coordinates')
        txs, local = annotate_contig(model, backbone, tokenizer, species, seqid, sequence, ordinal, out, 'cuda')
        transcripts.extend(txs); scopes.extend(local); lengths[seqid] = len(sequence)
    if not lengths:
        raise ValueError('Raw FASTA must contain at least one nonempty contig')
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen A43 emission model changed or received gradients')
    write_gff3(out/'A43.gff3', transcripts, lengths, 'M28_A43')
    D.R.C.R.write_summary(out, dict(arm='A43', completed=True, inference_complete=True,
        TRAIN_configuration=D.CONFIG, protocol=RAW_PROTOCOL, TRAIN_protocol=D.PROTOCOL,
        source_TRAIN_jobs=dict(fit=D.FIT_JOB, potentials=POTENTIAL_JOB, decode=DECODE_JOB, readout=READOUT_JOB, CPU=CHECK_JOB),
        checkpoint=fitted['final_checkpoint'], optimizer_steps_source=4608, model_parameters=366191,
        feature_revision=NATIVE.F.REVISION, feature_spec=D.R.C.R.FEATURE_SPEC,
        input_mode='explicit_raw_FASTA', input_fasta=str(fasta), lengths=lengths, scopes=scopes,
        context_policy=CONTEXT_POLICY, score_kind=SCORE_KIND, max_exons=128, threshold=0., null=0.,
        reference_used=False, labels_or_feature_cache_used=False, generator_model_loaded=False,
        query_or_beam_decoding=False, local_selected_only_owner_assembly=False,
        optimizer_steps_this_run=0, new_parameters=0, pretrained_encoder_updated=False,
        primary_CDS_chain_predictions=len(transcripts), wall_seconds=time.perf_counter()-began,
        limitations='Different base-context/whole-contig score and span domain from TRAIN windows; '
            'chromosome-scale resources, accuracy and independent generalization require their own evidence.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(), args.fasta.resolve())
