#!/usr/bin/env python3
"""A31 complete caller ruler; no reference before raw record/GFF reconciliation."""
import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import json
import math
from pathlib import Path
import infer_local_path_caller as I
import evaluate_chain_quality as A6
import evaluate_query_context as H
from evaluate_pilot import aggregate, stage_audit, validated_scope_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

A, E = A6.A, A6.A.E


def validate_complete(decoded, generator):
    expected = {(s, q, t): [[a, min(a+I.M.WIDTH, n)] for a in I.M.starts(n)]
                for s, qs in I.C.I.DEV_LENGTHS.items() for q, n in qs.items() for t in ('+', '-')}
    actual = {(s['species'], s['seqid'], s['strand']): s['windows'] for s in decoded['scopes']}
    if (decoded['arm'] != 'A31' or decoded['variant'] != 'A31_'+generator or
            decoded['completed'] is not True or decoded['inference_complete'] is not True or
            decoded['configuration'] != I.R.CONFIG or decoded['protocol'] != I.DEV_PROTOCOL or
            decoded['TRAIN_protocol'] != I.R.PROTOCOL or decoded['TRAIN_preparation_gate_passed'] is not True or
            decoded['generator_checkpoint'] != I.R.CONFIG['generators'][generator] or
            decoded['generator_optimizer_steps'] != 4608 or decoded['new_parameters'] != 0 or
            decoded['feature_revision'] != I.C.I.F.REVISION or decoded['feature_spec'] != I.C.R.FEATURE_SPEC or
            decoded['input_mode'] != 'fixed_DEV_raw_FASTA' or decoded['reference_used'] is not False or
            decoded['labels_or_feature_cache_used'] is not False or decoded['optimizer_steps_this_run'] != 0 or
            decoded['DEV_used_for_selection'] is not False or decoded['test_setaria_used'] is not False or
            decoded['score_kind'] != I.SCORE_KIND or decoded['windows'] != 8690 or
            len(decoded['scopes']) != 4 or actual != expected or
            decoded['lengths'] != {q: n for qs in I.C.I.DEV_LENGTHS.values() for q, n in qs.items()}):
        raise ValueError('Complete fixed 8690-window A31 raw-FASTA DEV grid required')


def paired_records(path, scope, old_path, old_scope, generator, counts):
    for row, old in zip_longest(validated_scope_records(path, scope), A6.records(old_path, old_scope)):
        if row is None or old is None:
            raise ValueError('Partial paired raw input scope')
        chains, scores, sources = row['chains'], row['score'], row['decoder_sources']
        if (row['variant'] != 'A31_'+generator or row['score_kind'] != I.SCORE_KIND or
                row['score_threshold'] != 0. or row['reference_used'] is not False or
                len(chains) != len(scores) or len(chains) != len(sources) or
                any(not math.isfinite(v) for v in scores) or
                any(not ss or not set(ss) <= {'greedy', 'beam4'} for ss in sources) or
                len({tuple(map(tuple, c)) for c in chains}) != len(chains)):
            raise ValueError('A31 natural pool, score or operating point differs')
        if any(row[k] != old[k] for k in ('window_index', 'window_id')):
            raise ValueError('Raw input window identity differs')
        if generator == 'A3' and any(row[k] != old[k] for k in ('chains', 'decoder_sources')):
            raise ValueError('Same frozen A3 natural geometry/order/source differs')
        # A23 has a different generator: do not demand or claim the same pool.
        actual = dict(free_candidates=len(chains), nonpositive_candidates=sum(v <= 0 for v in scores))
        if row['counts'] != actual:
            raise ValueError('A31 actual candidate/sign counts differ')
        counts.update(actual)
        yield row


def reconstruct(decoded, baseline, infer, generator):
    old_scopes = {(s['species'], s['seqid'], s['strand']): s for s in baseline['scopes']}
    predictions, assemblies, counts_by_scope = defaultdict(list), {}, {}
    for scope in decoded['scopes']:
        species, seqid, strand = scope['species'], scope['seqid'], scope['strand']
        old = old_scopes[(species, seqid, strand)]; counts = Counter()
        assembled = assemble_predictions(paired_records(infer/scope['file'], scope,
            I.C.INFER/old['file'], old, generator, counts), scope['windows'], strand)
        if dict(counts) != scope['candidate_counts'] or assembled['counts'] != scope['assembly_counts']:
            raise ValueError('A31 actual counts or assembly differ from inference')
        for index, row in enumerate(assembled['predictions']):
            ident = 'M28_A31_'+generator+'_'+str(scope['ordinal'])+'_'+strand+'_'+str(index)
            predictions[species].append(as_transcript(row['chain'], seqid, strand, ident, row['score']))
        key = species+'|'+strand; assemblies[key] = assembled; counts_by_scope[key] = dict(counts)
    transcripts = [t for values in predictions.values() for t in values]
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(infer/('A31_'+generator+'.gff3'), decoded['lengths']),
                                     complete_only=False)
    A.validate_gff_roundtrip(parsed, transcripts)
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final A31 GFF prediction count differs')
    return predictions, assemblies, counts_by_scope


def measure(predictions, assemblies, counts, context, infer, generator):
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
        raise ValueError('A31 selected-stage TP differs from exact-chain score')
    current = dict(per_species=per_species, pooled=pooled, stages=stages, candidate_counts=counts,
                   GFF3=str((infer/('A31_'+generator+'.gff3')).relative_to(I.ROOT)))
    H.check_ruler(current, current)  # New result must itself retain all 7728 references.
    return current


def main(infer, out, generator):
    prepared = I.R.completed(I.R.READOUT); checked = I.R.completed(I.R.CHECK)
    I.validate_train_preparation(prepared, checked, generator)
    if checked['source_jobs'] != {'readout': (I.R.READOUT/'JOBID').read_text().strip()}:
        raise ValueError('Actual source readout job must match CPU reconciliation')
    if not (I.ROOT/I.DEV_PROTOCOL).is_file():
        raise ValueError('Separate fixed A31 DEV protocol required')
    if out.exists():
        raise FileExistsError('A31 evaluation requires a new output directory')
    decoded = I.R.completed(infer); validate_complete(decoded, generator)
    old_fit = I.R.completed(I.C.FIT)
    baseline = I.R.completed(I.C.INFER); A6.validate_complete(old_fit, baseline)
    # Close actual pool, selection, geometry and GFF before reading any reference.
    predictions, assemblies, counts = reconstruct(decoded, baseline, infer, generator)
    current = measure(predictions, assemblies, counts, A.reference_context(), infer, generator)
    a3 = H.read_report(I.ROOT/'outputs/M28-A3-EVALUATION/summary.json')
    H.A3.Q.check_last_identity(a3['fit'])
    previous = H.read_report(I.ROOT/'outputs/M28-A-EVALUATION/summary.json')['methods']['A']
    control = H.read_report(I.ROOT/'reports/M28-METHOD-RESTART/r5_result.json')['methods']['C0']
    matched = H.read_report(I.ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json')
    if matched['completed'] is not True or matched['fit']['final_checkpoint'] != old_fit['final_checkpoint']:
        raise ValueError('Completed same A6 comparator required')
    native = H.read_report(I.ROOT/'outputs/M28-NATIVE-BASELINE-EVAL/summary.json')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Completed native comparison must retain unsupported reference strata')
    comparisons = dict(A3=a3['methods']['A3'], A6=matched['methods']['A6'], A=previous, R5_C0=control)
    for method in ('ANNEVO', 'Tiberius', 'Helixer'):
        comparisons[method] = dict(per_species={s: r['methods'][method] for s, r in native['species'].items()},
                                  pooled=native['pooled'][method])
    for values in comparisons.values():
        H.check_ruler(values, current)
    arm = 'A31_'+generator
    report = dict(experiment=arm, completed=True, protocol=I.DEV_PROTOCOL,
        status='complete_one_seed_bounded_DEV_not_independent_generalization',
        reference_primary_chains=7728, test_setaria_used=False,
        TRAIN_preparation=prepared, TRAIN_result_reconciliation=checked,
        inference_receipt={k: v for k, v in decoded.items() if k != 'scopes'}, methods={arm: current},
        same_raw_input_grid_and_native_four_block_context_verified=True,
        matched_A6_natural_candidate_geometry_order_sources_unchanged=True if generator == 'A3' else None,
        original_pilot_gate_auxiliary_only=A6.continuation_gate(current['per_species'], previous['per_species'], control['per_species']),
        contrasts={name: H.contrast(current, v) for name, v in comparisons.items()},
        exact_chain_comparison={name: dict(per_species={s: r['exact_chain'] for s, r in v['per_species'].items()},
            pooled=v['pooled']['exact_chain'], macro_F1=sum(r['exact_chain']['f1'] for r in v['per_species'].values())/2)
            for name, v in {arm: current, **comparisons}.items()},
        comparison_sources=dict(A3='outputs/M28-A3-EVALUATION/summary.json',
            A6='outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json', A='outputs/M28-A-EVALUATION/summary.json',
            R5_C0='reports/M28-METHOD-RESTART/r5_result.json', native='outputs/M28-NATIVE-BASELINE-EVAL/summary.json'),
        native_exposure_limit=native['source_limitations'],
        limitations='Repeated-use DEV only. A3 contrast changes pool and score; A6 has additional head exposure; '
            'A23 changes generator parameters and fitted shared weights. No matched-capacity causal claim, '
            'fair runtime comparison with cached A3, independent generalization or publication PASS.')
    out.mkdir(parents=True, exist_ok=False)
    I.C.R.write_summary(out, report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generator', choices=('A3', 'A23'), required=True)
    parser.add_argument('--inference-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); main(args.inference_dir.resolve(), args.output_dir.resolve(), args.generator)
