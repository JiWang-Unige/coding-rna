#!/usr/bin/env python3
"""A44 paired raw interface preparation; no actual raw execution is released."""
import argparse
import json
from pathlib import Path
import re
import time
import torch

import readout_local_structured_emission_train as D
import infer_natural_competitor as I
from src.m28.contig_emission_lattice import context_segments
from src.m28.c0_decode import genomic_chain
from src.m28.export import as_transcript, write_gff3
from src.m28.training import one_hot

ROOT, NATIVE, M, ARMS = D.ROOT, I.NATIVE, I.M, D.ARMS
POTENTIAL_JOB, DECODE_JOB, READOUT_JOB, CHECK_JOB = '13426150', '13426151', '13426152', '13426153'
RAW_PROTOCOL = 'reports/M28-METHOD-RESTART/a44_paired_raw_protocol.md'
CONTEXT_POLICY = I.CONTEXT_POLICY
SCORE_KIND = 'A44_paired_unit_role_gain_under_joint_base_owned_contig_context'
CPU_PROOFS = (
    'actual_both_original_draw_weights_loss_units_masks_and_4608_updates_verified',
    'actual_six_finite_paired_checkpoints_and_AdamW_exposure_verified',
    'actual_paired_last_model_potential_identity_and_complete_64_tasks_verified',
    'actual_new_DNA_ORF_motif_cap_span_and_full_window_FP64_log_probability_scores_verified',
    'actual_both_full_grammar_objectives_not_below_same_model_original_natural_pool',
    'actual_frozen_R4_unknowns_stronger_baseline_and_complete_denominators_recomputed',
    'actual_per_species_counts_and_fixed_three_five_condition_comparisons_recomputed',
)


def source_jobs():
    return dict(fits=D.FIT_JOBS, potentials=POTENTIAL_JOB, decode=DECODE_JOB,
                readout=READOUT_JOB, CPU=CHECK_JOB)


def validate_preparation(fitted, prepared, checked):
    if set(fitted) != set(ARMS):
        raise ValueError('Both matching A44 paired fits required')
    for arm in ARMS: D.validate_fit(fitted[arm], arm)
    source = dict(fits=D.FIT_JOBS, potentials=POTENTIAL_JOB, decode=DECODE_JOB)
    for report in (prepared, checked):
        if (report['completed'] is not True or report['configuration'] != D.CONFIG
                or report['protocol'] != D.PROTOCOL or report['windows_per_arm'] != 1491
                or report['original_draws'] != 1536 or report['original_candidates'] != 16321
                or report['reference_injected_predictions'] != 0 or report['optimizer_steps_this_run'] != 0
                or report['new_geometry_labels_computed'] is not True
                or any(report[k] is not False for k in ('model_or_backbone_forward', 'new_candidate_search',
                    'new_backbone_forward', 'R4_policy_or_catalog_changed', 'DEV_used',
                    'test_setaria_used', 'independent_species_used', 'A33_outputs_opened'))):
            raise ValueError('Actual complete original TRAIN-only paired readout and CPU required')
    if (prepared['source_jobs'] != source or checked['source_jobs'] != dict(source, readout=READOUT_JOB)
            or any(checked[k] is not True for k in CPU_PROOFS)
            or prepared['all_both_arm_natural_search_and_score_files_closed_before_labels'] is not True
            or prepared['original_R4_labels_exact'] is not True or checked['candidate_window_pairs'] != 19189):
        raise ValueError('Actual 64-task paired sources and seven CPU proofs required')
    counts, refs = prepared['selected_counts_per_species'], prepared['eligible_reference_pairs']
    expected = {'arabidopsis_thaliana': (1531, 347, 315, 548, 2172, 854),
                'oryza_sativa': (478, 251, 132, 152, 744, 232)}
    if (set(counts) != set(expected) or set(refs) != set(expected)
            or checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs):
        raise ValueError('Complete paired counts and original species denominators required')
    for species, (correct, nonmatch, unknown, long, eligible, five) in expected.items():
        baseline = dict(selected=correct+nonmatch+unknown, selected_correct=correct,
            selected_known_nonmatch=nonmatch, selected_unknown=unknown, selected_five_or_more_correct=long)
        if (counts[species]['A31_A3'] != baseline
                or refs[species] != dict(eligible=eligible, five_or_more=five)
                or set(counts[species]) != {'A31_A3', *ARMS}):
            raise ValueError('Actual stronger baseline and per-species full denominators required')
    gates, supported = D.gates(counts)
    if (prepared['predeclared_TRAIN_gates'] != gates or checked['predeclared_TRAIN_gates'] != gates
            or prepared['TRAIN_preparation_supported'] is not True
            or checked['TRAIN_preparation_supported'] is not True or not supported):
        raise ValueError('JOINT must actually pass both fixed five-condition comparisons')


def readiness():
    D.F.A.B.D.B.F.A.require_jobs(list(D.FIT_JOBS.values())+[POTENTIAL_JOB, READOUT_JOB, CHECK_JOB]
                               +[DECODE_JOB+'_'+str(i) for i in range(64)])
    fitted = {a: D.R.completed(D.F.output_root(a), D.FIT_JOBS[a]) for a in ARMS}
    prepared = D.R.completed(D.OUT, READOUT_JOB); checked = D.R.completed(D.CHECK, CHECK_JOB)
    validate_preparation(fitted, prepared, checked)
    return fitted


def load_models(fitted, device):
    models = {}
    for arm in ARMS:
        state = torch.load(ROOT/fitted[arm]['final_checkpoint'], map_location='cpu', weights_only=True)
        D.validate_fit(state, arm, state=True)
        model = D.F.A.A36.new_model(device).eval().requires_grad_(False)
        model.load_state_dict(state['model'], strict=True)
        if sum(p.numel() for p in model.parameters()) != 366191:
            raise ValueError('Only the two matching A44 single emission models allowed')
        models[arm] = model
    return models


@torch.inference_mode()
def paired_contig_emissions(models, backbone, tokenizer, sequence, locations, strand, device):
    """One native encoding per window, two models, original per-base ownership.

    Copy the same slices as stitch_contig_emissions without retaining every
    overlapping window or averaging roles. Minus row order changes, not columns.
    This prepares two independent outputs; it is not a model ensemble.
    """
    if strand not in ('+', '-') or set(sequence)-set('ACGTN') or set(models) != set(ARMS):
        raise ValueError('Normalized DNA, one strand and both separate A44 models required')
    segments = context_segments(len(sequence), locations)
    potentials = {a: torch.empty(len(sequence), 15, dtype=torch.float32) for a in ARMS}
    for index, (left, right) in enumerate(locations):
        raw = sequence[left:right]; oriented = raw if strand == '+' else M.rc(raw)
        features = NATIVE.encode_window(backbone, tokenizer, oriented, device)
        dna = one_hot(oriented, M.WIDTH)[None].to(device)
        _, start, end = segments[index]
        for arm in ARMS:
            logits = models[arm](features, dna, [right-left])['segmentation_logits'][0][:right-left].cpu()
            if (logits.dtype != torch.float32 or tuple(logits.shape) != (right-left, 15)
                    or not bool(torch.isfinite(logits).all())):
                raise ValueError('Finite valid-length paired raw fifteen-role FP32 logits required')
            if strand == '+': potentials[arm][start:end] = logits[start-left:end-left]
            else: potentials[arm][len(sequence)-end:len(sequence)-start] = logits[right-end:right-start]
        del features, dna, logits
    return dict(sequence=sequence if strand == '+' else M.rc(sequence), potentials=potentials,
                genomic_context_segments=segments)


@torch.inference_mode()
def annotate_contig(models, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    if not sequence: raise ValueError('Empty contig')
    locations = [[a, min(a+M.WIDTH, len(sequence))] for a in M.starts(len(sequence))]
    transcripts, scopes = {a: [] for a in ARMS}, {a: [] for a in ARMS}
    for strand in ('+', '-'):
        began = time.perf_counter()
        joined = paired_contig_emissions(models, backbone, tokenizer, sequence, locations, strand, device)
        encoding_seconds = time.perf_counter()-began
        for arm in ARMS:
            started = time.perf_counter()
            natural = D.D.A42.P.decode_emission_lattice(joined['sequence'], joined['potentials'][arm], D.CONFIG['max_exons'])
            natural['genomic_chains'] = tuple(tuple(genomic_chain(chain, 0, len(sequence), strand))
                                             for chain in natural['chains'])
            natural['genomic_context_segments'] = joined['genomic_context_segments']
            natural['context_policy'] = CONTEXT_POLICY
            path = out/arm/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand == '+' else 'm')+'.json')
            with path.open('x') as handle:
                json.dump(dict(variant=arm, score_kind=SCORE_KIND, species=species, seqid=seqid, strand=strand,
                    length=len(sequence), windows=locations, source_fit_job=D.FIT_JOBS[arm],
                    checkpoint=D.CONFIG['source_checkpoints'][arm], reference_used_in_scoring=False, **natural),
                    handle, allow_nan=False); handle.write('\n')
            for index, (chain, gain) in enumerate(zip(natural['genomic_chains'], natural['chain_gains'])):
                ident = 'M28_'+arm+'_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
                transcripts[arm].append(as_transcript(chain, seqid, strand, ident, gain))
            scopes[arm].append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
                windows=locations, file=str(path.relative_to(out)), arm=arm, context_policy=CONTEXT_POLICY,
                primary_CDS_chain_predictions=len(natural['chains']), decode_seconds=time.perf_counter()-started,
                shared_paired_encoding_seconds=encoding_seconds))
        del joined
    return transcripts, scopes


@torch.inference_mode()
def main(out, fasta):
    fitted = readiness()
    if not (ROOT/RAW_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed paired raw scientific/resource protocol before native execution')
    if fasta is None: raise ValueError('Explicit authorized FASTA required; no automatic DEV access')
    if out.exists(): raise FileExistsError('Paired raw inference requires a new output directory')
    models = load_models(fitted, 'cuda'); backbone, tokenizer = NATIVE.load_native()
    backbone.eval().requires_grad_(False)
    before = {a: {k: v.detach().clone() for k, v in m.state_dict().items()} for a, m in models.items()}
    out.mkdir(parents=True, exist_ok=False)
    for arm in ARMS: (out/arm).mkdir()
    began = time.perf_counter(); transcripts = {a: [] for a in ARMS}; scopes = {a: [] for a in ARMS}; lengths = {}
    for ordinal, (species, seqid, sequence) in enumerate(NATIVE.input_contigs(fasta)):
        if seqid in lengths: raise ValueError('Raw FASTA contig identifiers must be unique for GFF coordinates')
        txs, local = annotate_contig(models, backbone, tokenizer, species, seqid, sequence, ordinal, out, 'cuda')
        for arm in ARMS: transcripts[arm].extend(txs[arm]); scopes[arm].extend(local[arm])
        lengths[seqid] = len(sequence)
    if not lengths: raise ValueError('Raw FASTA must contain at least one nonempty contig')
    for arm, model in models.items():
        if (any(not torch.equal(v, before[arm][k]) for k, v in model.state_dict().items())
                or any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen paired raw model changed or received gradients')
        write_gff3(out/arm/(arm+'.gff3'), transcripts[arm], lengths, 'M28_'+arm)
    D.R.C.R.write_summary(out, dict(arms=list(ARMS), completed=True, inference_complete=True,
        TRAIN_configuration=D.CONFIG, protocol=RAW_PROTOCOL, TRAIN_protocol=D.PROTOCOL,
        source_TRAIN_jobs=source_jobs(), checkpoints=D.CONFIG['source_checkpoints'],
        optimizer_steps_source_per_arm=4608, model_parameters_per_arm=366191,
        feature_revision=NATIVE.F.REVISION, feature_spec=D.R.C.R.FEATURE_SPEC,
        input_mode='explicit_raw_FASTA', input_fasta=str(fasta), lengths=lengths, scopes=scopes,
        context_policy=CONTEXT_POLICY, score_kind=SCORE_KIND, max_exons=128, threshold=0., null=0.,
        reference_used=False, labels_or_feature_cache_used=False, generator_model_loaded=False,
        query_or_beam_decoding=False, local_selected_only_owner_assembly=False, models_ensembled=False,
        one_native_encoding_and_same_GEN_DNA_objects_per_paired_window=True,
        frozen_source_weights_unchanged=True, all_source_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0, new_parameters=0, pretrained_encoder_updated=False,
        primary_CDS_chain_predictions={a: len(t) for a, t in transcripts.items()}, wall_seconds=time.perf_counter()-began,
        limitations='Different base-context/whole-contig score and span domain from TRAIN windows; '
            'raw resources, common-ruler accuracy, independent generalization and novelty need their own evidence.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, required=True); parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(), args.fasta.resolve())
