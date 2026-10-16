#!/usr/bin/env python3
"""Three closed A3 DEV cells; saved-score diagnostic, never model selection."""
from collections import Counter
import csv
from itertools import zip_longest
import json
import math
import os
import subprocess
import evaluate_local_path_caller as V
import diagnose_local_caller_selector_bound as B
from evaluate_pilot import validated_scope_records
from src.m28.assembly import assemble_predictions

ROOT = B.ROOT
OUT = ROOT/'outputs/M28-A31-SCORE-POOL-DEV-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_score_pool_dev_diagnostic_protocol.md'
SOURCE_JOBS = dict(greedy_infer='13357568', greedy_eval='13357569',
                   union_infer='13399434', union_eval='13399454', references='13360116')
ARMS = ('greedy_query', 'greedy_local', 'union_local')


def paired_window(greedy, union):
    if any(greedy[k] != union[k] for k in ('window_index', 'window_id')):
        raise ValueError('Paired saved window identity differs')
    keys = [tuple(map(tuple, chain)) for chain in union['chains']]
    if (len(keys) != len(set(keys)) or len(keys) != len(union['score']) or
            len(keys) != len(union['decoder_sources'])):
        raise ValueError('Unique union geometry with one saved score/source per chain required')
    expected = {key for key, sources in zip(keys, union['decoder_sources']) if 'greedy' in sources}
    old = [tuple(map(tuple, chain)) for chain in greedy['chains']]
    if set(old) != expected:
        raise ValueError('Complete saved greedy geometry must equal the union greedy-source geometry')
    scores = dict(zip(keys, union['score']))
    if any(not math.isfinite(float(v)) for v in (*greedy['score'], *union['score'])):
        raise ValueError('Finite saved query and local scores required')
    # Keep original per-query duplicates; the unchanged assembly takes the max.
    return dict(window_index=greedy['window_index'], window_id=greedy['window_id'],
                chains=greedy['chains'], score=[scores[key] for key in old], reference_used=False)


def source_roots():
    jobs = list(SOURCE_JOBS.values())
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {}
    for line in result.stdout.splitlines():
        values = [s.strip() for s in line.split('|')]
        if len(values) >= 3 and values[0] in jobs:
            actual[values[0]] = values[1:3]
    if any(actual.get(job) != ['COMPLETED', '0:0'] for job in jobs):
        raise ValueError('Every original source root must close before output access')


def read_greedy(folder, job):
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or
            (folder/'JOBID').read_text().strip() != job):
        raise ValueError('Actual closed original A3 source identity required')
    return json.loads((folder/'summary.json').read_text())


def load_sources():
    source_roots()
    greedy_folder = ROOT/'outputs/M28-A3-DEV-INFER'
    old = read_greedy(greedy_folder, SOURCE_JOBS['greedy_infer'])
    old_eval = read_greedy(ROOT/'outputs/M28-A3-EVALUATION', SOURCE_JOBS['greedy_eval'])
    V.H.A3.validate_complete(old_eval['fit'], old)
    folder = ROOT/'outputs/M28-A31-A3-LOCAL-PATH-DEV-INFER'
    union = V.I.R.completed(folder, SOURCE_JOBS['union_infer'])
    measured = V.I.R.completed(ROOT/'outputs/M28-A31-A3-LOCAL-PATH-DEV-EVAL', SOURCE_JOBS['union_eval'])
    V.validate_complete(union, 'A3')
    if (old['checkpoint'] != union['generator_checkpoint'] or old['feature_revision'] != union['feature_revision'] or
            old_eval['reference_primary_chains'] != 7728 or old_eval['test_setaria_used'] or
            measured['experiment'] != 'A31_A3' or measured['protocol'] != V.I.DEV_PROTOCOL or
            measured['reference_primary_chains'] != 7728 or measured['test_setaria_used'] or
            measured['matched_A6_natural_candidate_geometry_order_sources_unchanged'] is not True):
        raise ValueError('Frozen generator, native qualification and complete original ruler required')
    return greedy_folder, old, old_eval['methods']['A3'], folder, union, measured['methods']['A31_A3']


def metrics(tp, predicted, reference):
    return dict(tp=tp, predicted=predicted, reference=reference,
                precision=tp/predicted if predicted else 0., recall=tp/reference,
                f1=2*tp/(predicted+reference))


def scope_result(assembled, references):
    keys = {r['chain'] for r in references}
    hits = {stage: len(set(chains) & keys) for stage, chains in assembled['stages'].items()}
    chosen = set(assembled['stages']['selected_unique'])
    return dict(reference=len(keys), exact_reference_hits=hits,
        not_in_actual_free_candidates=len(keys)-hits['free_unique'],
        lost_at_owner=hits['free_unique']-hits['owner_unique'],
        lost_at_conflict_or_nonpositive_gain=hits['owner_unique']-hits['selected_unique'],
        candidate_counts=assembled['counts'],
        selected_five_or_more_correct=sum(len(r['chain']) >= 5 and r['chain'] in chosen for r in references))


def contrast(old, new):
    return dict(count_changes={k: new[k]-old[k] for k in ('tp', 'predicted', 'selected_five_or_more_correct')},
                metric_change_percentage_points={k: 100*(new[k]-old[k]) for k in ('precision', 'recall', 'f1')})


def main():
    greedy_folder, old, old_measured, folder, union, union_measured = load_sources()
    OUT.mkdir(exist_ok=False)
    (OUT/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (OUT/'STATUS').write_text('RUNNING\n')
    old_scopes = {(s['species'], s['seqid'], s['strand']): s for s in old['scopes']}
    assembled_by_scope = {}; windows = 0
    with (OUT/'greedy_local_scores.jsonl').open('x') as output:
        for scope in union['scopes']:
            identity = (scope['species'], scope['seqid'], scope['strand'])
            greedy_scope = old_scopes[identity]
            if greedy_scope['windows'] != scope['windows']:
                raise ValueError('Paired full saved grid differs')
            counts = Counter()
            greedy = list(V.H.A3.query_scope_records(greedy_folder/greedy_scope['file'], greedy_scope, counts))
            V.H.A3.validate_diversity_counts(counts, greedy_scope)
            natural = list(validated_scope_records(folder/scope['file'], scope))
            middle = []
            for row, free in zip_longest(greedy, natural):
                if row is None or free is None:
                    raise ValueError('Partial paired saved population')
                if (free['variant'] != 'A31_A3' or free['score_kind'] != V.I.SCORE_KIND or
                        free['score_threshold'] != 0. or
                        any(not ss or not set(ss) <= {'greedy', 'beam4'} for ss in free['decoder_sources'])):
                    raise ValueError('Original fixed native scores and sources required')
                local = paired_window(row, free); middle.append(local); windows += 1
                output.write(json.dumps(local, allow_nan=False)+'\n')
            assembled = {arm: assemble_predictions(records, scope['windows'], scope['strand'])
                         for arm, records in zip(ARMS, (greedy, middle, natural))}
            if assembled['union_local']['counts'] != scope['assembly_counts']:
                raise ValueError('Original A31 union assembly must reproduce exactly')
            if assembled['greedy_query']['counts'] != old_measured['stages'][scope['species']+'|'+scope['strand']]['candidate_counts']:
                raise ValueError('Original A3 greedy assembly must reproduce exactly')
            for stage in ('free_unique', 'owner_unique'):
                if assembled['greedy_query']['stages'][stage] != assembled['greedy_local']['stages'][stage]:
                    raise ValueError('Scoring-only contrast must retain free and owner geometry')
            assembled_by_scope[identity] = assembled
    if windows != 8690:
        raise ValueError('All original 8690 paired saved windows required')
    # Geometry and reference-free assembly close before opening saved labels.
    saved = B.source(B.SOURCE, SOURCE_JOBS['references'])
    with (B.SOURCE/'references.tsv').open() as handle:
        references = B.prepare(list(csv.DictReader(handle, delimiter='\t')), saved)
    stages = {arm: {} for arm in ARMS}
    totals = {s: {a: Counter() for a in ARMS} for s in B.REFERENCES}
    with (OUT/'diagnostic_selections.jsonl').open('x') as output:
        for identity, assembled in assembled_by_scope.items():
            species, seqid, strand = identity
            local = [r for r in references if (r['species'], r['seqid'], r['strand']) == identity]
            for arm in ARMS:
                result = scope_result(assembled[arm], local)
                stages[arm][species+'|'+strand] = result
                totals[species][arm].update(tp=result['exact_reference_hits']['selected_unique'],
                    predicted=assembled[arm]['counts']['selected_unique'], reference=len(local),
                    selected_five_or_more_correct=result['selected_five_or_more_correct'])
                if arm != 'greedy_local':
                    original = (old_measured if arm == 'greedy_query' else union_measured)['stages'][species+'|'+strand]
                    if any(result[k] != original[k] for k in ('reference', 'exact_reference_hits',
                            'not_in_actual_free_candidates', 'lost_at_owner',
                            'lost_at_conflict_or_nonpositive_gain', 'candidate_counts')):
                        raise ValueError('Original complete source stage counts must reproduce')
            output.write(json.dumps(dict(species=species, seqid=seqid, strand=strand,
                greedy_local=assembled['greedy_local']['predictions']), allow_nan=False)+'\n')
    per_species = {s: {a: dict(metrics(c['tp'], c['predicted'], c['reference']),
        selected_five_or_more_correct=c['selected_five_or_more_correct']) for a, c in values.items()}
        for s, values in totals.items()}
    for species, total in B.REFERENCES.items():
        for arm in ARMS:
            if per_species[species][arm]['reference'] != total:
                raise ValueError('Retain all primary references including unsupported strata')
        for arm, measured in (('greedy_query', old_measured), ('union_local', union_measured)):
            if any(per_species[species][arm][k] != measured['per_species'][species]['exact_chain'][k]
                   for k in ('tp', 'predicted', 'reference')):
                raise ValueError('Both original cells must reproduce full source integer metrics')
    pooled = {}
    for arm in ARMS:
        sums = {k: sum(v[arm][k] for v in per_species.values())
                for k in ('tp', 'predicted', 'reference', 'selected_five_or_more_correct')}
        pooled[arm] = dict(metrics(sums['tp'], sums['predicted'], sums['reference']),
                           selected_five_or_more_correct=sums['selected_five_or_more_correct'])
    effects = {s: dict(score_change_within_greedy=contrast(v['greedy_query'], v['greedy_local']),
                      pool_expansion_at_local_score=contrast(v['greedy_local'], v['union_local']))
               for s, v in {**per_species, 'pooled': pooled}.items()}
    report = dict(completed=True, experiment='M28-A31-SCORE-POOL-DEV-DIAGNOSTIC', protocol=PROTOCOL,
        source_jobs=SOURCE_JOBS, windows=8690, reference_primary_chains=7728,
        generator_checkpoint=union['generator_checkpoint'], per_species=per_species, pooled=pooled,
        stages=stages, ordered_conditional_changes=effects,
        all_window_greedy_source_geometry_matched_before_reference_open=True,
        original_two_cells_full_stage_counts_and_integer_metrics_reproduced=True,
        original_greedy_duplicate_query_records_preserved=True, threshold=0.,
        new_model_forward=False, new_backbone_forward=False, new_candidate_search=False,
        new_R4_labelling=False, feature_cache_used=False, optimizer_steps_this_run=0,
        new_parameters=0, threshold_scan=False, test_setaria_used=False,
        diagnostic_intermediate_selection_recomputed=True, production_annotation_or_GFF_written=False,
        diagnostic_used_for_model_selection=False, original_TRAIN_gates_unchanged=True,
        limitation='Repeated-use DEV only. Ordered conditional contrasts, not causal shares or a factorial '
            'ablation. Cached original A3 versus native A31 geometry matches do not prove identical hidden '
            'activations. No union-query cell, interaction estimate, independent generalization or new gate.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
