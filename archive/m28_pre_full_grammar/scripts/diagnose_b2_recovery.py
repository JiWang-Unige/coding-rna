#!/usr/bin/env python3
"""Read-only B2 reference-chain recovery diagnosis after complete evaluation.

Outputs one TSV row per frozen primary reference and paired stage transitions.
No new inference, annotation parsing, filtering, tuning or mechanism gate.
Run real DEV processing through Slurm, after the main B2 terminal reconciliation.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
STAGES = ('free_unique', 'owner_unique', 'selected_unique')
ARMS = ('early', 'late')


def require_complete(decoded, evaluated):
    paired = evaluated.get('paired_inference_cost', {})
    if (evaluated.get('experiment') != 'M28-B2' or
        set(evaluated.get('methods', {})) != {'B2_early', 'B2_late'} or
        not decoded.get('inference_complete') or decoded.get('windows') != 8690 or
        decoded.get('reference_used') is not False or
        not decoded.get('control_replay', {}).get('discrete_exact') or
        not decoded.get('paired_common_gain', {}).get('exact_equality') or
        not decoded.get('paired_common_gain', {}).get('single_union_readout') or
        any(decoded.get(k) != paired.get(k)
            for k in ('checkpoint', 'base_checkpoint', 'feature_revision', 'windows'))):
        raise ValueError('Completed evaluation of this exact paired inference is required')


def chain_key(chain):
    return tuple(tuple(part[:2]) for part in chain)


def recovery_state(free, owner, selected, owner_gain):
    if selected:
        return 'selected'
    if not free:
        return 'not_free'
    if not owner:
        return 'lost_at_owner'
    if owner_gain is None:
        raise ValueError('Owned reference lacks its recorded gain')
    return 'owner_nonpositive' if owner_gain <= 0 else 'owner_positive_not_selected'


def reference_row(ref, stage_sets, owner_gains):
    chain = chain_key(ref['CDS'])
    row = {k: ref[k] for k in ('species', 'id', 'gene_id', 'seqid', 'strand')}
    row.update(exon_count=len(chain), span_bp=chain[-1][1]-chain[0][0],
               coding_bp=sum(b-a for a,b in chain),
               CDS_intervals=json.dumps(chain, separators=(',', ':')),
               unsupported_strata=';'.join(ref.get('M28_unsupported_strata_not_excluded', [])))
    for arm in ARMS:
        flags = [chain in stage_sets[arm][stage] for stage in STAGES]
        gain = owner_gains[arm].get(chain)
        for stage, flag in zip(STAGES, flags):
            row[arm+'_'+stage] = flag
        row[arm+'_owner_gain'] = gain
        row[arm+'_state'] = recovery_state(*flags, gain)
    return row


def transition_summary(records):
    """Counts are primary reference chains, never windows or all predictions."""
    result = {'reference_chains': len(records), 'stages': {}}
    for stage in STAGES:
        counts = Counter()
        for row in records:
            early, late = row['early_'+stage], row['late_'+stage]
            name = ('both' if early and late else 'early_only' if early else
                    'late_only' if late else 'neither')
            counts[name] += 1
        result['stages'][stage] = {name: counts[name]
                                  for name in ('both', 'early_only', 'late_only', 'neither')}
    states = Counter((r['late_state'], r['early_state']) for r in records)
    result['late_to_early_states'] = [
        {'late': late, 'early': early, 'reference_chains': count}
        for (late, early), count in sorted(states.items())]
    return result


def main(inference, evaluation, out):
    inference, evaluation, out = inference.resolve(), evaluation.resolve(), out.resolve()
    # The evaluator writes summary before its final artifact-cap check. A summary
    # alone is therefore insufficient evidence of successful terminal processing.
    if any((path/'STATUS').read_text().strip() != 'COMPLETED'
           for path in (inference, evaluation)):
        raise ValueError('Both paired inference and evaluation must be COMPLETED')
    decoded = json.loads((inference/'summary.json').read_text())
    evaluated = json.loads((evaluation/'summary.json').read_text())
    require_complete(decoded, evaluated)
    if out.exists():
        raise FileExistsError('Recovery diagnosis must use a new output directory')

    # Keep all real-data/dependency imports below the completion boundary.
    from evaluate_pilot import rows, validated_scope_records
    from src.m28.assembly import assemble_predictions, nearest_owner
    from src.m28.c0_decode import genomic_chain
    from fit_b2 import check_output_cap

    scopes = {arm: {(r['species'], r['seqid'], r['strand']): r
                    for r in decoded['arms'][arm]['scopes']} for arm in ARMS}
    expected = {(r['species'], r['seqid'], r['strand'])
                for r in evaluated['methods']['B2_early']['stages'].values()}
    if any(set(scopes[arm]) != expected for arm in ARMS):
        raise ValueError('Diagnostic scopes differ from the completed evaluation')
    species_names = sorted({s for s,_,_ in expected})
    references = {s: rows(ROOT/'outputs/M28-COMMON-RULER-R3/result'/(s+'.references.jsonl'))
                  for s in species_names}
    records = []
    for species, seqid, strand in sorted(expected):
        refs = [r for r in references[species] if r['seqid'] == seqid and r['strand'] == strand]
        reference = {chain_key(r['CDS']) for r in refs}
        stage_sets, gains = {}, {}
        for arm in ARMS:
            scope = scopes[arm][species, seqid, strand]
            windows = scope['windows']
            centers = [a+b for a,b in windows]
            gains[arm] = {}
            def observed_records():
                for record in validated_scope_records(inference/scope['file'], scope):
                    i = record['window_index']
                    a,b = windows[i]
                    for chain,score in zip(record['chains'], record['score']):
                        key = tuple(genomic_chain(chain, a, b, strand))
                        if key in reference and nearest_owner(key, windows, centers) == i:
                            gains[arm][key] = max(float(score), gains[arm].get(key, float('-inf')))
                    yield record
            assembled = assemble_predictions(observed_records(), windows, strand)
            stage_sets[arm] = {k: {chain_key(c) for c in assembled['stages'][k]} & reference for k in STAGES}
            counted = {k: len(stage_sets[arm][k]) for k in STAGES}
            official = evaluated['methods']['B2_'+arm]['stages'][species+'|'+strand]
            if counted != official['exact_reference_hits'] or assembled['counts'] != official['candidate_counts']:
                raise ValueError('Reference recovery replay differs from the official stage accounting')
            del assembled
        records.extend(reference_row(r, stage_sets, gains) for r in refs)

    if len(records) != evaluated['reference_primary_chains'] or len(records) != 7728:
        raise ValueError('Recovery diagnosis must retain the full primary reference denominator')
    summary = {
        'experiment': 'M28-B2-REFERENCE-RECOVERY',
        'scope': 'post-evaluation descriptive diagnosis; not a new gate or independent test',
        'inference_source': str(inference), 'evaluation_source': str(evaluation),
        'checkpoint': decoded['checkpoint'], 'test_setaria_used': False,
        'reference_filters_added': False, 'new_model_inference': False,
        'predeclared_gate_unchanged': evaluated['predeclared_mechanism_gate'],
        'state_definition': {
            'not_free': 'Absent from actual free candidates',
            'lost_at_owner': 'Free somewhere, absent from the assigned owner window',
            'owner_nonpositive': 'Owned, not selected, recorded owner gain <= 0',
            'owner_positive_not_selected': 'Owned with positive gain but not selected; no unique competitor cause claimed',
            'selected': 'Selected exact primary CDS chain'},
        'pooled': transition_summary(records),
        'per_species': {s: transition_summary([r for r in records if r['species'] == s])
                        for s in species_names}}

    out.mkdir()
    with (out/'references.tsv').open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    with (out/'summary.json').open('x') as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write('\n')
    if check_output_cap()+sum(p.stat().st_size for p in out.iterdir() if p.is_file()) > 2*1024**3:
        raise RuntimeError('Formal B2 plus diagnostic artifacts exceed the existing combined 2 GiB cap')
    print(json.dumps(summary['pooled'], indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inference-dir', type=Path, required=True)
    p.add_argument('--evaluation-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    args = p.parse_args()
    main(args.inference_dir, args.evaluation_dir, args.output_dir)
