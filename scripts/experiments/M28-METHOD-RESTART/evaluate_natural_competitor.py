"""A43 evaluation preparation; helpers only, no native/DEV execution entrypoint.

The eventual separately fixed execution contract must supply contig identities,
reference access and comparison receipts. Synthetic callers exercise this same
adapter without reading research inputs. Whole-contig outputs have no finite
proposal-pool/owner coverage statistic and contain no saved emission tensors.
"""
from collections import Counter
import json
import math
import re

import infer_natural_competitor as I
import evaluate_local_path_caller as V31
from evaluate_pilot import aggregate
from evaluate_structured_emission import comparison_table
from src.m28.c0_decode import genomic_chain
from src.m28.contig_emission_lattice import context_segments
from src.m28.export import as_transcript

A, E, A6, H = V31.A, V31.E, V31.A6, V31.H


def validate_complete(decoded, contigs):
    """Explicit seqid -> {species, length} mapping, not inferred from FASTA IDs."""
    expected_sources = dict(fit=I.D.FIT_JOB, potentials=I.POTENTIAL_JOB,
        decode=I.DECODE_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB)
    if (decoded['arm'] != 'A43' or decoded['completed'] is not True or
            decoded['inference_complete'] is not True or decoded['TRAIN_configuration'] != I.D.CONFIG or
            decoded['protocol'] != I.RAW_PROTOCOL or decoded['TRAIN_protocol'] != I.D.PROTOCOL or
            decoded['source_TRAIN_jobs'] != expected_sources or
            decoded['checkpoint'] != I.D.CONFIG['source_checkpoint'] or
            decoded['optimizer_steps_source'] != 4608 or decoded['model_parameters'] != 366191 or
            decoded['feature_revision'] != I.NATIVE.F.REVISION or decoded['feature_spec'] != I.D.R.C.R.FEATURE_SPEC or
            decoded['input_mode'] != 'explicit_raw_FASTA' or not decoded['input_fasta'] or
            decoded['context_policy'] != I.CONTEXT_POLICY or decoded['score_kind'] != I.SCORE_KIND or
            decoded['max_exons'] != 128 or decoded['threshold'] != 0. or decoded['null'] != 0. or
            decoded['optimizer_steps_this_run'] != 0 or decoded['new_parameters'] != 0 or
            any(decoded[k] is not False for k in ('reference_used', 'labels_or_feature_cache_used',
                'generator_model_loaded', 'query_or_beam_decoding', 'local_selected_only_owner_assembly',
                'pretrained_encoder_updated'))):
        raise ValueError('Complete matching single-model A43 raw inference receipt required')
    if not contigs or decoded['lengths'] != {q: c['length'] for q, c in contigs.items()}:
        raise ValueError('Explicit complete contig/species mapping and lengths required')
    scopes = decoded['scopes']
    expected = {(q, t) for q in contigs for t in ('+', '-')}
    if len(scopes) != len(expected) or {(s['seqid'], s['strand']) for s in scopes} != expected:
        raise ValueError('Complete nonduplicated both-strand contig scopes required')
    ordinals, files = {}, set()
    for s in scopes:
        q, t = s['seqid'], s['strand']; length = contigs[q]['length']
        grid = [[a, min(a+I.M.WIDTH, length)] for a in I.M.starts(length)]
        tag = 'p' if t == '+' else 'm'
        if (length <= 0 or not contigs[q]['species'] or s['species'] != 'user_input' or
                s['windows'] != grid or s['context_policy'] != I.CONTEXT_POLICY or
                type(s['ordinal']) is not int or s['ordinal'] < 0 or
                s['file'] != 'contig_'+str(s['ordinal']).zfill(4)+'_'+tag+'.json'):
            raise ValueError('Actual explicit-FASTA identity, original grid or joint context differs')
        ordinals.setdefault(q, s['ordinal'])
        if ordinals[q] != s['ordinal'] or s['file'] in files:
            raise ValueError('Saved contig ordinal/file identity differs')
        files.add(s['file'])
    if set(ordinals.values()) != set(range(len(contigs))):
        raise ValueError('Complete unique input contig ordinals required')


def validate_row(row, scope, length):
    chains, gains = row['chains'], row['chain_gains']
    if (row['variant'] != 'A43' or row['score_kind'] != I.SCORE_KIND or
            any(row[k] != scope[k] for k in ('species', 'seqid', 'strand', 'windows')) or
            row['length'] != length or row['reference_used'] is not False or
            row['reference_used_in_scoring'] is not False or row['context_policy'] != I.CONTEXT_POLICY or
            row['genomic_context_segments'] != [list(s) for s in context_segments(length, scope['windows'])] or
            len(chains) != len(gains) or len(chains) != scope['primary_CDS_chain_predictions'] or
            any(not math.isfinite(g) or g <= 0 for g in gains) or
            not math.isfinite(row['gain']) or
            (row['gain'] <= 0 if chains else row['gain'] != 0.) or
            not math.isclose(row['gain'], sum(gains), rel_tol=1e-10, abs_tol=1e-8)):
        raise ValueError('Actual joint-contig identity, positive chain gains, total or context differs')
    last_end = 0
    for chain in chains:
        if (not 1 <= len(chain) <= 128 or
                any(len(p) != 2 or any(type(v) is not int for v in p) or
                    not 0 <= p[0] < p[1] <= length for p in chain) or
                any(right[0]-left[1] < 5 for left, right in zip(chain, chain[1:])) or
                chain[0][0] < last_end):
            raise ValueError('Legal ordered nonoverlapping oriented chains, introns and exon cap required')
        last_end = chain[-1][1]
    converted = [[list(p) for p in genomic_chain(c, 0, length, scope['strand'])] for c in chains]
    if row['genomic_chains'] != converted:
        raise ValueError('Whole-contig oriented/genomic conversion or paired chain order differs')
    return converted, gains


def reconstruct(decoded, infer, contigs):
    """Close every saved chain and GFF coordinate/phase before reference access."""
    validate_complete(decoded, contigs)
    predictions = {c['species']: [] for c in contigs.values()}
    counts = {}
    for scope in decoded['scopes']:
        q, strand = scope['seqid'], scope['strand']; c = contigs[q]
        row = json.loads((infer/scope['file']).read_text())
        chains, gains = validate_row(row, scope, c['length'])
        for index, (chain, gain) in enumerate(zip(chains, gains)):
            ident = 'M28_A43_'+str(scope['ordinal'])+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
            predictions[c['species']].append(as_transcript(chain, q, strand, ident, gain))
        counts[q+'|'+strand] = dict(species=c['species'], seqid=q, strand=strand,
            selected_chains=len(chains), selected_gain=row['gain'])
    transcripts = [t for ts in predictions.values() for t in ts]
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final A43 prediction count differs')
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(infer/'A43.gff3', decoded['lengths']),
        complete_only=False)
    A.validate_gff_roundtrip(parsed, transcripts)
    return predictions, counts


def require_common_scope(contigs):
    if contigs != {q: dict(species=s, length=n) for s, q, n, _ in E.S.SCOPE}:
        raise ValueError('Complete original common-ruler contig/species mapping required')


def measure(predictions, counts, context, infer, contigs, input_sequences):
    """Unchanged full-reference/background ruler; never filter unreachable truth."""
    require_common_scope(contigs)
    refs, isoforms, parents, annotations, sequences, lengths = context
    if (lengths != {q: c['length'] for q, c in contigs.items()} or
            set(input_sequences) != set(contigs) or set(predictions) != set(A.REF_COUNTS) or
            set(refs) != set(A.REF_COUNTS)):
        raise ValueError('Complete original common-ruler species and contig context required')
    per_species = {}
    for species, seqid, length, _ in E.S.SCOPE:
        if (len({E.S.key(t) for t in refs[species]}) != A.REF_COUNTS[species] or
                len(sequences[species][seqid]) != length):
            raise ValueError('Full reference denominator or raw sequence length differs; do not filter truth')
        dna = re.sub('[^ACGT]', 'N', sequences[species][seqid].upper())
        if re.sub('[^ACGT]', 'N', input_sequences[seqid].upper()) != dna:
            raise ValueError('Reported raw FASTA differs from the complete reference contig DNA')
        failures = Counter()
        for tx in predictions[species]:
            failures.update(I.M.chain_checks(tx, dna)[0])
        if failures:
            raise ValueError('A43 exported predictions violate unchanged DNA/ORF/phase/splice grammar: '+str(dict(failures)))
        measured = A6.score_species(predictions[species], refs[species], isoforms[species],
            length, annotations[species], parents[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        measured['CDS_structure_diagnostic'] = dict(fully_supported=len(predictions[species]),
            predictions=len(predictions[species]), failure_counts_nonexclusive={})
        per_species[species] = measured
    current = dict(per_species=per_species, pooled=aggregate(per_species),
        selected_scope_counts=counts, GFF3=str(infer/'A43.gff3'))
    H.check_ruler(current, current)
    return current


def compare(current, previous, alternative, native):
    for old in (previous, alternative, native):
        if old['reference_primary_chains'] != 7728 or old['test_setaria_used'] is not False:
            raise ValueError('Completed original full-ruler comparison receipts required')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Native comparisons must retain unsupported reference strata')
    methods = dict(A31_A3=previous['methods']['A31_A3'], A31_A23=alternative['methods']['A31_A23'])
    for name in ('ANNEVO', 'Tiberius', 'Helixer'):
        methods[name] = dict(per_species={s: r['methods'][name] for s, r in native['species'].items()},
            pooled=native['pooled'][name])
    return dict(exact_chain_comparison=comparison_table(dict(A43=current, **methods), current),
        contrasts={name: H.contrast(current, v) for name, v in methods.items()},
        native_exposure_limit=native['source_limitations'])


def evaluate_saved(decoded, infer, contigs, input_loader, reference_loader, comparison_loader):
    """Caller owns the eventual fixed execution scope; this module reads no DEV by itself."""
    require_common_scope(contigs)
    predictions, counts = reconstruct(decoded, infer, contigs)
    input_sequences = input_loader(decoded['input_fasta'])
    if (set(input_sequences) != set(contigs) or
            any(len(input_sequences[q]) != c['length'] for q, c in contigs.items())):
        raise ValueError('Complete explicitly mapped raw FASTA sequences required before reference access')
    current = measure(predictions, counts, reference_loader(), infer, contigs, input_sequences)
    compared = compare(current, *comparison_loader())
    return dict(completed=True, methods=dict(A43=current), **compared,
        reference_primary_chains=7728, test_setaria_used=False,
        raw_predictions_and_GFF_reconciled_before_reference_access=True,
        saved_chain_coordinates_phases_positive_gains_and_total_reconciled=True,
        reported_raw_FASTA_matches_complete_reference_contig_DNA=True,
        actual_DNA_ORF_phase_splice_and_exon_cap_verified=True,
        unsupported_reference_strata_retained=True, finite_pool_owner_coverage_applicable=False,
        independent_neural_forward_replay=False, independent_saved_emission_score_replay=False,
        independent_full_grammar_DP_replay=False,
        status='complete_one_seed_repeated_DEV_not_independent_generalization',
        limitations='Caller must preserve the fixed FASTA between inference and evaluation. Saved outputs '
            'contain no emission tensors: gain totals are not independent neural/DP '
            'score replay. Different base context and span domain from TRAIN/old owner assembly; '
            'no matched-capacity causal claim, fair native runtime, independent generalization or publication PASS. '
            'Nonexact primary predictions are not automatically biological false positives.')
