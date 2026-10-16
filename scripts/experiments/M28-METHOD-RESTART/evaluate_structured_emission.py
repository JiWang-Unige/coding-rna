#!/usr/bin/env python3
"""A36 full DEV ruler, only after saved raw predictions/GFF reconcile."""
import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import math
from pathlib import Path
import infer_structured_emission as I
import evaluate_local_path_caller as V31
from evaluate_pilot import aggregate, stage_audit, validated_scope_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

A, E, A6, H = V31.A, V31.E, V31.A6, V31.H
BASELINE_INFER = I.ROOT/'outputs/M28-A31-A3-LOCAL-PATH-DEV-INFER'
BASELINE_EVAL = I.ROOT/'outputs/M28-A31-A3-LOCAL-PATH-DEV-EVAL'
A23_EVAL = I.ROOT/'outputs/M28-A31-A23-LOCAL-PATH-DEV-EVAL'


def validate_complete(decoded):
    expected = {(s, q, t): [[a, min(a+I.M.WIDTH, n)] for a in I.M.starts(n)]
                for s, qs in I.C.I.DEV_LENGTHS.items() for q, n in qs.items() for t in ('+', '-')}
    actual = {(s['species'], s['seqid'], s['strand']): s['windows'] for s in decoded['scopes']}
    if (decoded['arm'] != 'A36' or decoded['completed'] is not True or decoded['inference_complete'] is not True or
            decoded['configuration'] != I.D.CONFIG or decoded['protocol'] != I.DEV_PROTOCOL or
            decoded['TRAIN_protocol'] != I.D.PROTOCOL or decoded['fit_protocol'] != I.D.F.PROTOCOL or
            decoded['TRAIN_preparation_gate_passed'] is not True or
            decoded['source_TRAIN_jobs'] != dict(fit=I.D.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB) or
            decoded['generator_checkpoint'] != I.C.GENERATOR or
            decoded['emission_checkpoint'] != str((I.D.F.OUT/'step_004608.pt').relative_to(I.ROOT)) or
            decoded['generator_optimizer_steps'] != 4608 or decoded['emission_optimizer_steps'] != 4608 or
            decoded['model_parameters'] != dict(generator=649457, emission=366191) or
            decoded['feature_revision'] != I.C.I.F.REVISION or decoded['feature_spec'] != I.C.R.FEATURE_SPEC or
            decoded['input_mode'] != 'fixed_DEV_raw_FASTA' or decoded['reference_used'] is not False or
            decoded['labels_or_feature_cache_used'] is not False or decoded['optimizer_steps_this_run'] != 0 or
            decoded['DEV_used_for_selection'] is not False or decoded['test_setaria_used'] is not False or
            decoded['A33_outputs_opened'] is not False or
            decoded['GEN_encoded_once_per_oriented_window'] is not True or
            decoded['natural_proposal_generator'] != 'unchanged_A3_greedy_plus_beam4' or
            decoded['model_weights_unchanged'] is not True or
            decoded['all_parameters_frozen_and_gradients_absent'] is not True or
            decoded['score_kind'] != I.SCORE_KIND or decoded['windows'] != 8690 or
            len(decoded['scopes']) != 4 or actual != expected or
            decoded['lengths'] != {q: n for qs in I.C.I.DEV_LENGTHS.values() for q, n in qs.items()}):
        raise ValueError('Complete fixed 8690-window two-model A36 raw-FASTA DEV grid required')


def paired_records(path, scope, old_path, old_scope, counts):
    for row, old in zip_longest(validated_scope_records(path, scope), validated_scope_records(old_path, old_scope)):
        if row is None or old is None:
            raise ValueError('Partial paired raw input scope')
        chains, scores, sources = row['chains'], row['score'], row['decoder_sources']
        if (row['variant'] != 'A36' or row['score_kind'] != I.SCORE_KIND or
                row['score_threshold'] != 0. or row['reference_used'] is not False or
                len(chains) != len(scores) or len(chains) != len(sources) or
                any(not math.isfinite(v) for v in scores) or
                any(not ss or not set(ss) <= {'greedy', 'beam4'} for ss in sources) or
                len({tuple(map(tuple, c)) for c in chains}) != len(chains)):
            raise ValueError('A36 natural pool, score or operating point differs')
        if (old['variant'] != 'A31_A3' or old['score_kind'] != V31.I.SCORE_KIND or
                old['score_threshold'] != 0. or old['reference_used'] is not False):
            raise ValueError('Actual stronger A31_A3 raw candidate source required')
        if any(row[k] != old[k] for k in ('window_index', 'window_id', 'chains', 'decoder_sources')):
            raise ValueError('Same frozen A3 natural geometry/order/source differs')
        actual = dict(free_candidates=len(chains), nonpositive_candidates=sum(v <= 0 for v in scores))
        if row['counts'] != actual:
            raise ValueError('Actual A36 candidate/sign counts differ')
        counts.update(actual)
        yield row


def reconstruct(decoded, baseline, infer, baseline_infer):
    old_scopes = {(s['species'], s['seqid'], s['strand']): s for s in baseline['scopes']}
    predictions, assemblies, counts_by_scope = defaultdict(list), {}, {}
    for scope in decoded['scopes']:
        species, seqid, strand = scope['species'], scope['seqid'], scope['strand']
        old = old_scopes[(species, seqid, strand)]; counts = Counter()
        assembled = assemble_predictions(paired_records(infer/scope['file'], scope,
            baseline_infer/old['file'], old, counts), scope['windows'], strand)
        if dict(counts) != scope['candidate_counts'] or assembled['counts'] != scope['assembly_counts']:
            raise ValueError('Actual A36 candidate counts or saved assembly differ')
        for index, row in enumerate(assembled['predictions']):
            ident = 'M28_A36_'+str(scope['ordinal'])+'_'+('p' if strand == '+' else 'm')+'_'+str(index)
            predictions[species].append(as_transcript(row['chain'], seqid, strand, ident, row['score']))
        key = species+'|'+strand; assemblies[key] = assembled; counts_by_scope[key] = dict(counts)
    transcripts = [t for values in predictions.values() for t in values]
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(infer/'A36.gff3', decoded['lengths']), complete_only=False)
    A.validate_gff_roundtrip(parsed, transcripts)
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final A36 GFF prediction count differs')
    return predictions, assemblies, counts_by_scope


def measure(predictions, assemblies, counts, context, infer):
    refs, isoforms, parents, annotations, sequences, _ = context
    per_species, stages = {}, {}
    for key, assembled in assemblies.items():
        species, strand = key.split('|')
        audit = stage_audit(assembled, [r for r in refs[species] if r['strand'] == strand])
        audit.update(species=species, seqid=next(iter(I.C.I.DEV_LENGTHS[species])), strand=strand)
        stages[key] = audit
    for species, seqid, length, _ in E.S.SCOPE:
        measured = A6.score_species(predictions[species], refs[species], isoforms[species],
                                   length, annotations[species], parents[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        failures = Counter(); valid = 0
        for tx in predictions[species]:
            problems, _ = I.M.chain_checks(tx, sequences[species][seqid])
            failures.update(problems); valid += not problems
        measured['CDS_structure_diagnostic'] = dict(fully_supported=valid,
            predictions=len(predictions[species]), failure_counts_nonexclusive=dict(failures))
        per_species[species] = measured
    pooled = aggregate(per_species)
    pooled['actual_free_primary_coverage'] = sum(x['exact_reference_hits']['free_unique'] for x in stages.values())/7728
    pooled['owner_primary_coverage'] = sum(x['exact_reference_hits']['owner_unique'] for x in stages.values())/7728
    if sum(x['exact_reference_hits']['selected_unique'] for x in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('Selected-stage TP differs from exact-chain score')
    current = dict(per_species=per_species, pooled=pooled, stages=stages, candidate_counts=counts,
                   GFF3=str((infer/'A36.gff3').relative_to(I.ROOT)))
    H.check_ruler(current, current)
    return current


def comparison_table(methods, current):
    table = {}
    for name, values in methods.items():
        H.check_ruler(values, current)
        per_species = {s: r['exact_chain'] for s, r in values['per_species'].items()}
        table[name] = dict(per_species=per_species, pooled=values['pooled']['exact_chain'],
            macro_F1=sum(r['f1'] for r in per_species.values())/len(per_species),
            worst_species_F1=min(r['f1'] for r in per_species.values()))
    return table


def main(infer, out, source_job):
    I.readiness()
    if not (I.ROOT/I.DEV_PROTOCOL).is_file():
        raise ValueError('Separate fixed A36 DEV protocol required')
    if out.exists():
        raise FileExistsError('A36 evaluation requires a new output directory')
    I.D.F.A.require_jobs([str(source_job), '13399434', '13399454', '13399455'])
    decoded = I.D.R.completed(infer, source_job); validate_complete(decoded)
    baseline = I.D.R.completed(BASELINE_INFER, '13399434'); V31.validate_complete(baseline, 'A3')
    previous = I.D.R.completed(BASELINE_EVAL, '13399454')
    alternative = I.D.R.completed(A23_EVAL, '13399455')
    if (previous['reference_primary_chains'] != 7728 or previous['test_setaria_used'] or
            alternative['reference_primary_chains'] != 7728 or alternative['test_setaria_used'] or
            previous['inference_receipt'] != {k: v for k, v in baseline.items() if k != 'scopes'}):
        raise ValueError('Actual stronger baseline evaluation or complete ruler differs')
    # Rebuild all actual prediction records/GFF before opening reference inputs.
    predictions, assemblies, counts = reconstruct(decoded, baseline, infer, BASELINE_INFER)
    current = measure(predictions, assemblies, counts, A.reference_context(), infer)
    a3 = H.read_report(I.ROOT/'outputs/M28-A3-EVALUATION/summary.json')
    H.A3.Q.check_last_identity(a3['fit'])
    a6 = H.read_report(I.ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json')
    if a6['completed'] is not True:
        raise ValueError('Completed A6 comparator required')
    control = H.read_report(I.ROOT/'reports/M28-METHOD-RESTART/r5_result.json')['methods']['C0']
    historical = H.read_report(I.ROOT/'outputs/M28-A-EVALUATION/summary.json')['methods']['A']
    native = H.read_report(I.ROOT/'outputs/M28-NATIVE-BASELINE-EVAL/summary.json')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Completed native comparisons must retain unsupported reference strata')
    comparisons = dict(A31_A3=previous['methods']['A31_A3'], A31_A23=alternative['methods']['A31_A23'],
        A3=a3['methods']['A3'], A6=a6['methods']['A6'], A=historical, R5_C0=control)
    for method in ('ANNEVO', 'Tiberius', 'Helixer'):
        comparisons[method] = dict(per_species={s: r['methods'][method] for s, r in native['species'].items()},
                                  pooled=native['pooled'][method])
    table = comparison_table(dict(A36=current, **comparisons), current)
    report = dict(experiment='M28-A36-STRUCTURED-EMISSION-DEV', completed=True, protocol=I.DEV_PROTOCOL,
        status='complete_one_seed_repeated_DEV_not_independent_generalization',
        reference_primary_chains=7728, test_setaria_used=False, source_inference_job=str(source_job),
        source_TRAIN_jobs=decoded['source_TRAIN_jobs'], inference_receipt={k: v for k, v in decoded.items() if k != 'scopes'},
        methods=dict(A36=current), same_raw_input_grid_and_native_four_block_context_verified=True,
        matched_A31_A3_natural_candidate_geometry_order_sources_unchanged=True,
        raw_predictions_and_GFF_reconciled_before_reference_access=True,
        exact_chain_comparison=table, contrasts={name: H.contrast(current, v) for name, v in comparisons.items()},
        comparison_sources=dict(A31_A3=str((BASELINE_EVAL/'summary.json').relative_to(I.ROOT)),
            A31_A23=str((A23_EVAL/'summary.json').relative_to(I.ROOT)),
            A3='outputs/M28-A3-EVALUATION/summary.json', A6='outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json',
            A='outputs/M28-A-EVALUATION/summary.json', R5_C0='reports/M28-METHOD-RESTART/r5_result.json',
            native='outputs/M28-NATIVE-BASELINE-EVAL/summary.json'), native_exposure_limit=native['source_limitations'],
        limitation='Repeated two-species DEV only. A36 adds representation capacity and supervised exposure; '
            'A31_A3 shares proposals, A31_A23 changes generator. No matched-capacity causal claim, '
            'fair native runtime comparison with cached methods, independent generalization or publication PASS.')
    out.mkdir(parents=True, exist_ok=False); I.C.R.write_summary(out, report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inference-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--source-job', required=True)
    args = parser.parse_args(); main(args.inference_dir.resolve(), args.output_dir.resolve(), args.source_job)
