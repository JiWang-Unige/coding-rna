"""A44 paired whole-contig evaluator helpers; no raw/DEV execution entrypoint."""
import json
import math

import infer_local_structured_emission as I
import evaluate_natural_competitor as V
from src.m28.c0_decode import genomic_chain
from src.m28.contig_emission_lattice import context_segments
from src.m28.export import as_transcript


def validate_complete(decoded, contigs):
    if (decoded['arms'] != list(I.ARMS) or decoded['completed'] is not True
            or decoded['inference_complete'] is not True or decoded['TRAIN_configuration'] != I.D.CONFIG
            or decoded['protocol'] != I.RAW_PROTOCOL or decoded['TRAIN_protocol'] != I.D.PROTOCOL
            or decoded['source_TRAIN_jobs'] != I.source_jobs()
            or decoded['checkpoints'] != I.D.CONFIG['source_checkpoints']
            or decoded['optimizer_steps_source_per_arm'] != 4608 or decoded['model_parameters_per_arm'] != 366191
            or decoded['feature_revision'] != I.NATIVE.F.REVISION or decoded['feature_spec'] != I.D.R.C.R.FEATURE_SPEC
            or decoded['input_mode'] != 'explicit_raw_FASTA' or not decoded['input_fasta']
            or decoded['context_policy'] != I.CONTEXT_POLICY or decoded['score_kind'] != I.SCORE_KIND
            or decoded['max_exons'] != 128 or decoded['threshold'] != 0. or decoded['null'] != 0.
            or decoded['optimizer_steps_this_run'] != 0 or decoded['new_parameters'] != 0
            or any(decoded[k] is not True for k in ('one_native_encoding_and_same_GEN_DNA_objects_per_paired_window',
                'frozen_source_weights_unchanged', 'all_source_parameters_frozen_and_gradients_absent'))
            or any(decoded[k] is not False for k in ('reference_used', 'labels_or_feature_cache_used',
                'generator_model_loaded', 'query_or_beam_decoding', 'local_selected_only_owner_assembly',
                'models_ensembled', 'pretrained_encoder_updated'))):
        raise ValueError('Complete matching separate two-model A44 raw receipt required')
    if (not contigs or decoded['lengths'] != {q: c['length'] for q, c in contigs.items()}
            or set(decoded['scopes']) != set(I.ARMS)
            or set(decoded['primary_CDS_chain_predictions']) != set(I.ARMS)):
        raise ValueError('Explicit complete paired contig/species mapping and both arms required')
    expected = {(q, t) for q in contigs for t in ('+', '-')}; paired_ordinals = None
    for arm in I.ARMS:
        scopes = decoded['scopes'][arm]
        if len(scopes) != len(expected) or {(s['seqid'], s['strand']) for s in scopes} != expected:
            raise ValueError('Complete nonduplicated both-strand scopes for each arm required')
        ordinals, files = {}, set()
        for s in scopes:
            q, strand = s['seqid'], s['strand']; length = contigs[q]['length']
            grid = [[a, min(a+I.M.WIDTH, length)] for a in I.M.starts(length)]
            tag = 'p' if strand == '+' else 'm'
            if (length <= 0 or not contigs[q]['species'] or s['species'] != 'user_input' or s['arm'] != arm
                    or s['windows'] != grid or s['context_policy'] != I.CONTEXT_POLICY
                    or type(s['ordinal']) is not int or s['ordinal'] < 0
                    or s['file'] != arm+'/contig_'+str(s['ordinal']).zfill(4)+'_'+tag+'.json'):
                raise ValueError('Actual paired arm, explicit-FASTA identity, original grid or context differs')
            ordinals.setdefault(q, s['ordinal'])
            if ordinals[q] != s['ordinal'] or s['file'] in files:
                raise ValueError('Saved paired ordinal/file identity differs')
            files.add(s['file'])
        if set(ordinals.values()) != set(range(len(contigs))) or (paired_ordinals is not None and ordinals != paired_ordinals):
            raise ValueError('Complete identical input contig ordinals for both arms required')
        paired_ordinals = ordinals


def validate_row(row, scope, length, arm):
    chains, gains = row['chains'], row['chain_gains']
    if (row['variant'] != arm or row['source_fit_job'] != I.D.FIT_JOBS[arm]
            or row['checkpoint'] != I.D.CONFIG['source_checkpoints'][arm] or row['score_kind'] != I.SCORE_KIND
            or any(row[k] != scope[k] for k in ('species', 'seqid', 'strand', 'windows'))
            or row['length'] != length or row['reference_used'] is not False
            or row['reference_used_in_scoring'] is not False or row['context_policy'] != I.CONTEXT_POLICY
            or row['genomic_context_segments'] != [list(s) for s in context_segments(length, scope['windows'])]
            or len(chains) != len(gains) or len(chains) != scope['primary_CDS_chain_predictions']
            or any(not math.isfinite(g) or g <= 0 for g in gains)
            or not math.isfinite(row['gain']) or (row['gain'] <= 0 if chains else row['gain'] != 0.)
            or not math.isclose(row['gain'], sum(gains), rel_tol=1e-10, abs_tol=1e-8)):
        raise ValueError('Actual paired arm/source/context, positive chain gains or total differs')
    last_end = 0
    for chain in chains:
        if (not 1 <= len(chain) <= 128 or any(len(p) != 2 or any(type(v) is not int for v in p)
                or not 0 <= p[0] < p[1] <= length for p in chain)
                or any(right[0]-left[1] < 5 for left, right in zip(chain, chain[1:])) or chain[0][0] < last_end):
            raise ValueError('Legal ordered nonoverlapping oriented chains, introns and exon cap required')
        last_end = chain[-1][1]
    converted = [[list(p) for p in genomic_chain(c, 0, length, scope['strand'])] for c in chains]
    if row['genomic_chains'] != converted:
        raise ValueError('Paired whole-contig conversion or chain/gain order differs')
    return converted, gains


def reconstruct(decoded, infer, contigs):
    """Both complete arms and GFF files close before FASTA/reference loaders."""
    validate_complete(decoded, contigs); predictions, counts = {}, {}
    for arm in I.ARMS:
        predictions[arm] = {c['species']: [] for c in contigs.values()}; counts[arm] = {}
        for scope in decoded['scopes'][arm]:
            q, strand = scope['seqid'], scope['strand']; c = contigs[q]
            row = json.loads((infer/scope['file']).read_text()); chains, gains = validate_row(row, scope, c['length'], arm)
            for index, (chain, gain) in enumerate(zip(chains, gains)):
                ident = 'M28_'+arm+'_'+str(scope['ordinal'])+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
                predictions[arm][c['species']].append(as_transcript(chain, q, strand, ident, gain))
            counts[arm][q+'|'+strand] = dict(species=c['species'], seqid=q, strand=strand,
                selected_chains=len(chains), selected_gain=row['gain'])
        transcripts = [t for ts in predictions[arm].values() for t in ts]
        if len(transcripts) != decoded['primary_CDS_chain_predictions'][arm]:
            raise ValueError('Final paired prediction count differs')
        parsed = V.E.S.E.primary_transcripts(V.E.S.E.parse_annotation(infer/arm/(arm+'.gff3'), decoded['lengths']),
                                            complete_only=False)
        V.A.validate_gff_roundtrip(parsed, transcripts)
    return predictions, counts


def compare(currents, previous, alternative, native):
    for old in (previous, alternative, native):
        if old['reference_primary_chains'] != 7728 or old['test_setaria_used'] is not False:
            raise ValueError('Completed original full-ruler comparison receipts required')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Native comparisons must retain unsupported reference strata')
    methods = dict(A31_A3=previous['methods']['A31_A3'], A31_A23=alternative['methods']['A31_A23'])
    for name in ('ANNEVO', 'Tiberius', 'Helixer'):
        methods[name] = dict(per_species={s: r['methods'][name] for s, r in native['species'].items()},
            pooled=native['pooled'][name])
    return dict(exact_chain_comparison=V.comparison_table(dict(**currents, **methods), currents['A44_JOINT']),
        contrasts={a: {name: V.H.contrast(current, v) for name, v in methods.items()} for a, current in currents.items()},
        paired_JOINT_vs_LOCAL=V.H.contrast(currents['A44_JOINT'], currents['A44_LOCAL']),
        native_exposure_limit=native['source_limitations'])


def evaluate_saved(decoded, infer, contigs, input_loader, reference_loader, comparison_loader):
    """Helpers only: the eventual fixed caller owns access and execution scope."""
    V.require_common_scope(contigs)
    predictions, counts = reconstruct(decoded, infer, contigs)
    inputs = input_loader(decoded['input_fasta'])
    if set(inputs) != set(contigs) or any(len(inputs[q]) != c['length'] for q, c in contigs.items()):
        raise ValueError('Complete explicitly mapped raw FASTA required before reference access')
    context = reference_loader(); currents = {}
    for arm in I.ARMS:
        # Reuse unchanged full-primary/background arithmetic, not A43 identity.
        current = V.measure(predictions[arm], counts[arm], context, infer/arm, contigs, inputs)
        current['GFF3'] = str(infer/arm/(arm+'.gff3')); currents[arm] = current
    compared = compare(currents, *comparison_loader())
    return dict(completed=True, methods=currents, **compared, source_TRAIN_jobs=I.source_jobs(),
        reference_primary_chains=7728, test_setaria_used=False,
        both_raw_prediction_and_GFF_files_reconciled_before_reference_access=True,
        actual_complete_common_reference_and_background_used_for_both_arms=True,
        reported_raw_FASTA_matches_complete_reference_contig_DNA=True,
        actual_DNA_ORF_phase_splice_and_exon_cap_verified=True,
        unsupported_reference_strata_retained=True, finite_pool_owner_coverage_applicable=False,
        models_ensembled=False, independent_neural_forward_replay=False,
        independent_saved_emission_score_replay=False, independent_full_grammar_DP_replay=False,
        status='complete_one_seed_repeated_DEV_not_independent_generalization',
        limitations='Caller must preserve fixed FASTA and access scope. Two separate outputs share GEN, '
            'not an ensemble. Saved outputs contain no emissions for independent neural/DP replay; '
            'base-context/span domain differs from TRAIN windows. No matched total supervision, '
            'fair runtime, independent generalization, novelty or publication PASS; nonexact primary '
            'predictions are not automatically biological false positives.')
