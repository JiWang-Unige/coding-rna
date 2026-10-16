#!/usr/bin/env python3
"""Split closed A31 native DEV attrition; no new scores, predictions or tuning."""
from collections import Counter, defaultdict
import csv
import json
import math
import os
import subprocess
import diagnose_completed_quality_scores as D
import diagnose_local_caller_selector_bound as B
import evaluate_local_path_caller as V
from evaluate_pilot import validated_scope_records
from src.m28.assembly import nearest_owner
from src.m28.c0_decode import genomic_chain
from src.m28.core import select_nonoverlapping

ROOT = B.ROOT
OUT = ROOT/'outputs/M28-A31-DEV-SCORE-ATTRITION'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_dev_score_attrition_protocol.md'
SOURCE_JOBS = dict(A3_infer='13399434', A3_eval='13399454',
                   A23_infer='13399441', A23_eval='13399455', A6_attrition='13360116')
ARMS = ('A6', 'A31_A3', 'A31_A23')


def source_roots():
    jobs = list(SOURCE_JOBS.values())
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {}
    for line in result.stdout.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) >= 3 and parts[0] in jobs:
            actual[parts[0]] = parts[1:3]
    if any(actual.get(job) != ['COMPLETED', '0:0'] for job in jobs):
        raise ValueError('All actual source roots must close before opening A31 outputs')


def load_sources():
    source_roots()
    sources = {}
    for generator in ('A3', 'A23'):
        infer = ROOT/('outputs/M28-A31-'+generator+'-LOCAL-PATH-DEV-INFER')
        evaluation = ROOT/('outputs/M28-A31-'+generator+'-LOCAL-PATH-DEV-EVAL')
        decoded = V.I.R.completed(infer, SOURCE_JOBS[generator+'_infer'])
        measured = V.I.R.completed(evaluation, SOURCE_JOBS[generator+'_eval'])
        V.validate_complete(decoded, generator)
        if (measured['experiment'] != 'A31_'+generator or measured['protocol'] != V.I.DEV_PROTOCOL or
                measured['reference_primary_chains'] != 7728 or measured['test_setaria_used'] or
                not measured['same_raw_input_grid_and_native_four_block_context_verified'] or
                measured['matched_A6_natural_candidate_geometry_order_sources_unchanged'] is not
                    (True if generator == 'A3' else None)):
            raise ValueError('Closed original A31 native ruler and pool qualification required')
        sources[generator] = (infer, decoded, measured['methods']['A31_'+generator])
    saved = B.source(B.SOURCE, SOURCE_JOBS['A6_attrition'])
    with (B.SOURCE/'references.tsv').open() as handle:
        references = B.prepare(list(csv.DictReader(handle, delimiter='\t')), saved)
    return sources, references


def scope_pool(records, scope):
    """Replay the original owner/max-duplicate/positive-gain assembly rule."""
    windows = scope['windows']; centers = [a+b for a, b in windows]
    free, owned = set(), {}
    total = owner_records = 0
    for record in records:
        index = record['window_index']; a, b = windows[index]
        if len(record['chains']) != len(record['score']):
            raise ValueError('One saved score required per candidate')
        for chain, value in zip(record['chains'], record['score']):
            key = tuple(genomic_chain(chain, a, b, scope['strand'])); gain = float(value)
            if not math.isfinite(gain):
                raise ValueError('Nonfinite saved native gain')
            total += 1; free.add(key)
            if nearest_owner(key, windows, centers) == index:
                owner_records += 1
                owned[key] = max(owned.get(key, -math.inf), gain)
    keys = sorted(owned)
    selected = {keys[i] for i in select_nonoverlapping(keys, [owned[k] for k in keys])}
    actual = dict(free_records=total, free_unique=len(free), owner_records=owner_records,
                  owner_unique=len(owned), selected_unique=len(selected))
    if actual != {k: scope['assembly_counts'][k] for k in actual}:
        raise ValueError('Saved original native assembly counts must reproduce')
    return free, owned, selected


def scope_states(references, free, owned, selected, audit):
    labels = {r['id']: D.state(r['chain'], free, owned, r['chain'] in selected) for r in references}
    counts = Counter(labels.values())
    actual = (len(labels), counts['not_free'], counts['lost_owner'],
              counts['nonpositive_gain']+counts['positive_unselected'], counts['selected'])
    expected = (audit['reference'], audit['not_in_actual_free_candidates'], audit['lost_at_owner'],
                audit['lost_at_conflict_or_nonpositive_gain'], audit['exact_reference_hits']['selected_unique'])
    if actual != expected:
        raise ValueError('Every reference and original combined-stage count must reproduce')
    return labels


def main():
    sources, references = load_sources()
    if OUT.exists():
        raise FileExistsError('Preserve any earlier diagnostic output')
    diagnostics = {(r['species'], r['id']): dict(r) for r in references}
    counts = {s: {arm: Counter() for arm in ARMS} for s in B.REFERENCES}
    strata = {s: defaultdict(lambda: {arm: Counter() for arm in ARMS}) for s in B.REFERENCES}
    transitions = {s: {arm: Counter() for arm in ARMS[1:]} for s in B.REFERENCES}
    for generator, (folder, decoded, measured) in sources.items():
        arm = 'A31_'+generator
        for scope in decoded['scopes']:
            species, seqid, strand = scope['species'], scope['seqid'], scope['strand']
            local = [r for r in references if (r['species'], r['seqid'], r['strand']) ==
                     (species, seqid, strand)]
            free, owned, selected = scope_pool(validated_scope_records(folder/scope['file'], scope), scope)
            labels = scope_states(local, free, owned, selected, measured['stages'][species+'|'+strand])
            for row in local:
                key = (species, row['id']); label = labels[row['id']]
                if generator == 'A3' and (
                        (label == 'not_free') != (row['A6_state'] == 'not_free') or
                        (label == 'lost_owner') != (row['A6_state'] == 'lost_owner')):
                    raise ValueError('Matched A31 A3 and A6 free/owner membership must be identical')
                diagnostics[key][arm+'_state'] = label
                diagnostics[key][arm+'_owner_gain'] = owned.get(row['chain'])
                transitions[species][arm][row['A6_state']+'->'+label] += 1
    for row in diagnostics.values():
        species = row['species']; row['CDS_intervals'] = json.dumps(row.pop('chain'))
        size = len(json.loads(row['CDS_intervals']))
        group = 'one' if size == 1 else 'two_to_four' if size <= 4 else 'five_or_more'
        for arm in ARMS:
            label = row[arm+'_state']
            counts[species][arm][label] += 1; strata[species][group][arm][label] += 1
    normalized = {s: {a: {k: c[k] for k in D.STATES} for a, c in arms.items()} for s, arms in counts.items()}
    for species, total in B.REFERENCES.items():
        if any(sum(values.values()) != total for values in normalized[species].values()):
            raise ValueError('All unsupported and supported references must remain')
        for generator, (_, _, measured) in sources.items():
            if normalized[species]['A31_'+generator]['selected'] != measured['per_species'][species]['exact_chain']['tp']:
                raise ValueError('Reference selection must match actual native exact-chain TP')
    OUT.mkdir(exist_ok=False)
    (OUT/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (OUT/'STATUS').write_text('RUNNING\n')
    with (OUT/'references.tsv').open('x', newline='') as handle:
        fields = list(next(iter(diagnostics.values())))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        writer.writerows(diagnostics.values())
    report = dict(completed=True, experiment='M28-A31-DEV-SCORE-ATTRITION', protocol=PROTOCOL,
        source_jobs=SOURCE_JOBS, reference_primary_chains=7728, per_species=normalized,
        pooled={a: {k: sum(v[a][k] for v in normalized.values()) for k in D.STATES} for a in ARMS},
        per_species_exon_class={s: {g: {a: {k: c[k] for k in D.STATES} for a, c in arms.items()}
            for g, arms in groups.items()} for s, groups in strata.items()},
        reference_state_transitions_vs_A6=transitions,
        actual_saved_owner_scores_selection_full_denominators_and_original_stage_counts_reproduced=True,
        new_training=False, new_model_forward=False, new_backbone_forward=False,
        new_predictions=False, new_R4_labelling=False, reference_filters_added=False,
        optimizer_steps_this_run=0, threshold_scan=False, tau=0., test_setaria_used=False,
        interpretation='Repeated-use DEV diagnostic, not independent generalization or causal attribution. '
            'Positive unselected is an observed selector exclusion, not proof of mutually conflicting true genes. '
            'A31 A3 shares the A6 pool but not its head or training exposure; A23 changes generator and scores. '
            'No new gate, F1, threshold choice, annotation or release.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
