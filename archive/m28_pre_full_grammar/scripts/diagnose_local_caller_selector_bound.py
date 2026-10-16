#!/usr/bin/env python3
"""Label-assisted geometric ceilings from closed A6 pools, not new annotations."""
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import subprocess
from src.m28.core import select_nonoverlapping

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/'outputs/M28-A6-A7-SCORE-ATTRITION'
NATIVE = ROOT/'outputs/M28-NATIVE-BASELINE-EVAL'
OUT = ROOT/'outputs/M28-A31-FIXED-POOL-SELECTOR-BOUND'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_fixed_pool_selector_bound_protocol.md'
SOURCE_JOBS = dict(saved_attrition='13360116', native_evaluation='13358432')
REFERENCES = dict(arabidopsis_thaliana=5437, oryza_sativa=2291)
STATES = ('not_free', 'lost_owner', 'nonpositive_gain', 'positive_unselected', 'selected')
MODES = ('all_primary', 'A6_free', 'A6_owner')


def source(folder, job):
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or
            (folder/'JOBID').read_text().strip() != job):
        raise ValueError('Actual closed source identity required')
    value = json.loads((folder/'summary.json').read_text())
    if not value['completed'] or value['reference_primary_chains'] != 7728 or value['test_setaria_used']:
        raise ValueError('Complete unchanged 7728-reference source required')
    return value


def source_roots():
    jobs = list(SOURCE_JOBS.values())
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {}
    for line in result.stdout.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) >= 3 and parts[0] in jobs:
            actual[parts[0]] = parts[1:3]
    if any(actual.get(j) != ['COMPLETED', '0:0'] for j in jobs):
        raise ValueError('Successful actual source roots required')


def oracle(rows):
    """Unit gains maximize exact-reference cardinality under the unchanged selector."""
    scopes = defaultdict(list)
    for row in rows:
        scopes[(row['species'], row['seqid'], row['strand'])].append(row)
    selected = []
    for values in scopes.values():
        indices = select_nonoverlapping([r['chain'] for r in values], [1.]*len(values))
        selected.extend(values[i] for i in indices)
    return selected


def comparison(tp, reference, baseline):
    if baseline['reference'] != reference:
        raise ValueError('Same full reference denominator required for the native comparison')
    btp, bpred = baseline['tp'], baseline['predicted']
    delta = tp*(bpred+reference)-btp*(tp+reference)
    return dict(max_span_nonoverlap_TP=tp, full_reference_denominator=reference,
        label_assisted_precision_one_F1_upper=2*tp/(reference+tp),
        ideal_F1_strictly_above_native_ANNEVO=delta > 0,
        minimum_TP_at_precision_one_to_strictly_exceed_native_ANNEVO=
            (btp*reference)//(bpred+reference-btp)+1,
        maximum_nonexact_predictions_if_TP_equals_oracle_and_F1_strictly_exceeds_ANNEVO=
            (delta-1)//btp if delta > 0 else None)


def prepare(rows, saved):
    counts = {s: {a: Counter() for a in ('A6', 'A7')} for s in REFERENCES}
    identities, keys, result = set(), set(), []
    for row in rows:
        species = row['species']; chain = tuple(tuple(p) for p in json.loads(row['CDS_intervals']))
        ident = (species, row['id']); key = (species, row['seqid'], row['strand'], chain)
        if ident in identities or key in keys:
            raise ValueError('Full unique primary identities and exact-chain keys required')
        identities.add(ident); keys.add(key)
        for arm in ('A6', 'A7'):
            label = row[arm+'_state']
            if label not in STATES or row[arm+'_selected'] not in ('True', 'False') or \
                    (row[arm+'_selected'] == 'True') != (label == 'selected'):
                raise ValueError('Actual closed pool/selection membership differs')
            counts[species][arm][label] += 1
        result.append(dict(species=species, seqid=row['seqid'], strand=row['strand'],
            id=row['id'], chain=chain, A6_state=row['A6_state'], unsupported_strata=row['unsupported_strata']))
    normalized = {s: {a: {k: c[k] for k in STATES} for a, c in arms.items()} for s, arms in counts.items()}
    if normalized != saved['per_species'] or any(sum(v['A6'].values()) != REFERENCES[s]
                                                for s, v in normalized.items()):
        raise ValueError('All original source states and full species denominators must reproduce')
    return result


def main():
    source_roots()
    saved = source(SOURCE, SOURCE_JOBS['saved_attrition'])
    native = source(NATIVE, SOURCE_JOBS['native_evaluation'])
    if (saved['reference_filters_added'] or saved['new_training'] or saved['new_model_inference'] or
            saved['new_predictions'] or saved['optimizer_steps_this_run'] != 0 or
            not native['unsupported_reference_strata_retained']):
        raise ValueError('Accepted closed source-only labels and native ruler required')
    with (SOURCE/'references.tsv').open() as handle:
        rows = prepare(list(csv.DictReader(handle, delimiter='\t')), saved)
    populations = dict(all_primary=rows,
        A6_free=[r for r in rows if r['A6_state'] != 'not_free'],
        A6_owner=[r for r in rows if r['A6_state'] not in ('not_free', 'lost_owner')])
    measured, selected_ids = {}, {}
    for mode, population in populations.items():
        selected = oracle(population); chosen = {(r['species'], r['id']) for r in selected}
        selected_ids[mode] = chosen; per_species = {}
        for species, reference in REFERENCES.items():
            available = [r for r in population if r['species'] == species]
            tp = sum(r['species'] == species for r in selected)
            per_species[species] = dict(correct_candidates=len(available),
                correct_candidates_lost_only_to_span_conflict=len(available)-tp,
                **comparison(tp, reference, native['species'][species]['methods']['ANNEVO']['exact_chain']))
        tp = len(selected)
        measured[mode] = dict(per_species=per_species, pooled=dict(correct_candidates=len(population),
            correct_candidates_lost_only_to_span_conflict=len(population)-tp,
            **comparison(tp, 7728, native['pooled']['ANNEVO']['exact_chain'])))
    OUT.mkdir(exist_ok=False)
    (OUT/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (OUT/'STATUS').write_text('RUNNING\n')
    with (OUT/'references.tsv').open('x', newline='') as handle:
        fields = ('species', 'id', 'seqid', 'strand', 'CDS_intervals', 'unsupported_strata',
                  'A6_state', *[m+'_oracle_selected' for m in MODES])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for row in rows:
            writer.writerow(dict(species=row['species'], id=row['id'], seqid=row['seqid'], strand=row['strand'],
                CDS_intervals=json.dumps(row['chain']), unsupported_strata=row['unsupported_strata'],
                A6_state=row['A6_state'], **{m+'_oracle_selected': (row['species'], row['id']) in chosen
                                           for m, chosen in selected_ids.items()}))
    report = dict(completed=True, experiment='M28-A31-FIXED-POOL-SELECTOR-BOUND', protocol=PROTOCOL,
        source_jobs=SOURCE_JOBS, reference_primary_chains=7728, measured=measured,
        actual_closed_state_memberships_and_full_denominators_reproduced=True,
        oracle_unit_scores_use_reference_labels=True, reference_filters_added=False,
        new_model_forward=False, new_backbone_forward=False, optimizer_steps_this_run=0,
        current_A31_outputs_opened=False, new_annotation_or_GFF=False, test_setaria_used=False,
        limitation='Label-assisted geometric ceilings, not learned attainable scores or model predictions. '
            'A6 owner bound constrains A31_A3 only if its later actual same-pool check passes; '
            'not A31_A23/A30, independent causality, fair timing, background gate or generalization.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
