#!/usr/bin/env python3
"""Read-only A44 saved TRAIN gains, finite conflicts and natural retention."""
import argparse
from collections import Counter, defaultdict
import json
import math
import statistics
import time
import torch

import infer_local_structured_emission as I
from src.m28.phase_evidence_chain import candidate_roles
from src.m28.structured_emission import trainable_path_scores

D, ROOT, ARMS = I.D, I.ROOT, I.ARMS
OUT = ROOT/'outputs/M28-A44-PAIRED-MECHANISM'
CHECKS = ROOT/'outputs/M28-A44-PAIRED-MECHANISM-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a44_paired_mechanism_diagnostic.md'
ROLE_GROUPS = dict(coding_phase=(1, 2, 3), intron_phase=(4, 5, 6),
                   donor_phase=(7, 8, 9), acceptor_phase=(10, 11, 12), start=(13,), stop=(14,))


def key(chain):
    return tuple(tuple(exon) for exon in chain)


def overlap(a, b):
    return a[0][0] < b[-1][1] and b[0][0] < a[-1][1]


def validate_sources(prepared, checked):
    """Numeric closure is required; a scientific FAIL remains diagnosable."""
    source = dict(fits=D.FIT_JOBS, potentials=I.POTENTIAL_JOB, decode=I.DECODE_JOB)
    for value, jobs in ((prepared, source), (checked, dict(source, readout=I.READOUT_JOB))):
        if (value['completed'] is not True or value['source_jobs'] != jobs
                or value['configuration'] != D.CONFIG or value['protocol'] != D.PROTOCOL
                or value['windows_per_arm'] != 1491 or value['original_draws'] != 1536
                or value['original_candidates'] != 16321 or value['reference_injected_predictions'] != 0
                or value['new_geometry_labels_computed'] is not True
                or any(value[k] != v for k, v in D.SCOPE.items())):
            raise ValueError('Actual closed original A44 TRAIN sources required')
    if (checked['candidate_window_pairs'] != 19189
            or any(checked[k] is not True for k in I.CPU_PROOFS)
            or prepared['all_both_arm_natural_search_and_score_files_closed_before_labels'] is not True
            or prepared['original_R4_labels_exact'] is not True
            or any(prepared[k] != checked[k] for k in ('selected_counts_per_species',
                'eligible_reference_pairs', 'predeclared_TRAIN_gates', 'TRAIN_preparation_supported'))):
        raise ValueError('Both closed decisions and actual seven numeric proofs required')
    gates, supported = D.gates(prepared['selected_counts_per_species'])
    if gates != prepared['predeclared_TRAIN_gates'] or supported != prepared['TRAIN_preparation_supported']:
        raise ValueError('Preserve the actual frozen gates, including scientific FAIL')


def readiness(checks_job):
    jobs = list(D.FIT_JOBS.values())+[I.POTENTIAL_JOB, I.READOUT_JOB, I.CHECK_JOB, checks_job]
    jobs += [I.DECODE_JOB+'_'+str(i) for i in range(64)]
    D.F.A.B.D.B.F.A.require_jobs(jobs)
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED'
            or (CHECKS/'JOBID').read_text().strip() != checks_job):
        raise ValueError('Actual artificial mechanism checks required')
    prepared = D.R.completed(D.OUT, I.READOUT_JOB)
    checked = D.R.completed(D.CHECK, I.CHECK_JOB)
    validate_sources(prepared, checked)
    return prepared, checked


def components(local, chain):
    """Original FP64 role-minus-intergenic sums, not probabilities or a new score."""
    roles = candidate_roles(chain, len(local), local.device)
    start, end = chain[0][0], chain[-1][1]
    positions = torch.arange(start, end, device=local.device)
    chosen = roles[start:end]
    logits = local[positions].double()
    gains = logits[torch.arange(end-start, device=local.device), chosen]-logits[:, 0]
    totals = torch.bincount(chosen, weights=gains, minlength=15)
    values = {name: float(totals[list(indices)].sum()) for name, indices in ROLE_GROUPS.items()}
    if not math.isclose(sum(values.values()), float(gains.sum()), rel_tol=1e-12, abs_tol=1e-8):
        raise ValueError('Role decomposition must conserve the unchanged whole-chain gain')
    return values


@torch.inference_mode()
def diagnose_window(candidates, potentials, predictions):
    chains = [row['chain'] for row in candidates]
    scores = {a: trainable_path_scores(potentials[a], chains).tolist() for a in ARMS}
    selected = {a: {key(c) for c in predictions[a]['chains']} for a in ARMS}
    nonmatches = [i for i, row in enumerate(candidates) if row['label'] == 0]
    rows = []
    for i, candidate in enumerate(candidates):
        if candidate['label'] != 1:
            continue
        chain = candidate['chain']; conflicts = [j for j in nonmatches if overlap(chain, chains[j])]
        retained = {a: key(chain) in selected[a] for a in ARMS}
        state = ('both_selected' if all(retained.values()) else 'newly_selected' if retained[ARMS[1]]
                 else 'lost_selected' if retained[ARMS[0]] else 'neither_selected')
        arms = {}
        for arm in ARMS:
            best = max(conflicts, key=lambda j: (scores[arm][j], -j)) if conflicts else None
            terms = components(potentials[arm], chain)
            if not math.isclose(sum(terms.values()), scores[arm][i], rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError('Paired role terms differ from the original unit chain score')
            arms[arm] = dict(gain=scores[arm][i], gain_positive=scores[arm][i] > 0., selected=retained[arm],
                components=terms, finite_best_known_conflict_index=best,
                finite_best_known_conflict_gain=scores[arm][best] if best is not None else None,
                finite_known_conflict_margin=scores[arm][i]-scores[arm][best] if best is not None else None)
        old, new = (arms[a] for a in ARMS)
        rows.append(dict(candidate_index=i, chain=chain, exons=len(chain), five_or_more=len(chain) >= 5,
            selection_transition=state, finite_known_conflict_count=len(conflicts), arms=arms,
            joint_minus_local_gain=new['gain']-old['gain'],
            joint_minus_local_components={k: new['components'][k]-old['components'][k] for k in ROLE_GROUPS},
            joint_minus_local_finite_known_conflict_margin=(
                new['finite_known_conflict_margin']-old['finite_known_conflict_margin'] if conflicts else None)))
    return rows


def distribution(values):
    values = list(values)
    return dict(n=len(values), mean=statistics.fmean(values) if values else None,
                median=statistics.median(values) if values else None,
                minimum=min(values) if values else None, maximum=max(values) if values else None,
                positive=sum(v > 0 for v in values), negative=sum(v < 0 for v in values),
                zero=sum(v == 0 for v in values))


def summarize(rows):
    counts = Counter(row['selection_transition'] for row in rows)
    result = dict(reference_pairs=len(rows), five_or_more_reference_pairs=sum(r['five_or_more'] for r in rows),
        selection_transitions={s: counts[s] for s in ('both_selected', 'newly_selected', 'lost_selected', 'neither_selected')},
        references_without_finite_known_conflict=sum(not r['finite_known_conflict_count'] for r in rows),
        joint_minus_local_gain=distribution(r['joint_minus_local_gain'] for r in rows),
        joint_minus_local_finite_known_conflict_margin=distribution(
            r['joint_minus_local_finite_known_conflict_margin'] for r in rows
            if r['joint_minus_local_finite_known_conflict_margin'] is not None),
        joint_minus_local_components={k: distribution(r['joint_minus_local_components'][k] for r in rows)
                                     for k in ROLE_GROUPS}, arms={})
    for arm in ARMS:
        result['arms'][arm] = dict(selected=sum(r['arms'][arm]['selected'] for r in rows),
            selected_five_or_more=sum(r['arms'][arm]['selected'] and r['five_or_more'] for r in rows),
            positive_not_selected=sum(r['arms'][arm]['gain_positive'] and not r['arms'][arm]['selected'] for r in rows),
            nonpositive_not_selected=sum(not r['arms'][arm]['gain_positive'] and not r['arms'][arm]['selected'] for r in rows))
    return result


@torch.inference_mode()
def main(checks_job):
    prepared, checked = readiness(checks_job); began = time.perf_counter()
    potentials = D.require_potentials(I.POTENTIAL_JOB)
    _, _, support, _ = D.F.P.A.population('13424089')
    natural = {a: D.R.C.R.rows(D.OUT/(a+'_natural_scores.jsonl')) for a in ARMS}
    if any(len(natural[a]) != 1491 for a in ARMS):
        raise ValueError('Both complete original saved predictions required')
    groups = defaultdict(list); unknowns = Counter(); finite_known = Counter()
    with (OUT/'reference_pairs.jsonl').open('x') as handle:
        for ordinal in range(1491):
            values = {a: D.load_potential(a, ordinal, potentials['records'][a][ordinal], potentials) for a in ARMS}
            predictions = {a: natural[a][ordinal] for a in ARMS}
            w = values[ARMS[0]]['window']; source = support[w['id']]
            if (source['original_ordinal'] != ordinal or source['window'] != w
                    or any(values[a]['window'] != w or predictions[a]['window'] != w
                           or predictions[a]['arm'] != a or predictions[a]['original_ordinal'] != ordinal
                           for a in ARMS)):
                raise ValueError('Same original oriented TRAIN window and arm required')
            s = w['species']; candidates = source['candidates']
            unknowns[s] += sum(c['label'] == -1 for c in candidates)
            finite_known[s] += sum(c['label'] >= 0 for c in candidates)
            rows = diagnose_window(candidates, {a: values[a]['local_logits'] for a in ARMS}, predictions)
            for row in rows:
                row.update(window_id=w['id'], original_ordinal=ordinal, species=s, strand=w['strand'])
                groups[s].append(row); handle.write(json.dumps(row, allow_nan=False)+'\n')
            del values
            if (ordinal+1) % 128 == 0 or ordinal+1 == 1491:
                print(json.dumps(dict(saved_paired_windows_diagnosed=ordinal+1, total=1491)), flush=True)
    summaries = {s: summarize(rows) for s, rows in groups.items()}
    if set(summaries) != set(prepared['eligible_reference_pairs']):
        raise ValueError('Complete original species required')
    for s, value in summaries.items():
        refs = prepared['eligible_reference_pairs'][s]
        if (value['reference_pairs'] != refs['eligible'] or value['five_or_more_reference_pairs'] != refs['five_or_more']
                or any(value['arms'][a]['selected'] != prepared['selected_counts_per_species'][s][a]['selected_correct']
                       or value['arms'][a]['selected_five_or_more'] != prepared['selected_counts_per_species'][s][a]['selected_five_or_more_correct']
                       for a in ARMS)):
            raise ValueError('Paired exact gold identity must reproduce complete readout denominators and selections')
    if sum(len(rows) for rows in groups.values()) != 2916 or sum(unknowns.values())+sum(finite_known.values()) != 19189:
        raise ValueError('Original complete gold and fixed support population required')
    D.R.C.R.write_summary(OUT, dict(completed=True, protocol=PROTOCOL, source_configuration=D.CONFIG,
        source_jobs=I.source_jobs(), checks_job=checks_job, windows=1491, reference_pairs=2916,
        candidate_window_pairs=19189, finite_known_candidates_by_species=dict(finite_known),
        unknown_candidates_excluded_from_negative_conflicts_by_species=dict(unknowns),
        per_species=summaries, by_selection_transition={s: {state: summarize([r for r in rows if r['selection_transition'] == state])
            for state in ('both_selected', 'newly_selected', 'lost_selected', 'neither_selected')} for s, rows in groups.items()},
        original_predeclared_TRAIN_gates=prepared['predeclared_TRAIN_gates'],
        original_TRAIN_preparation_supported=checked['TRAIN_preparation_supported'],
        diagnostic_changes_gates=False, model_weights_opened=False, model_or_backbone_forward=False,
        new_candidate_search=False, new_geometry_labels_computed=False, seconds=time.perf_counter()-began,
        limitation='Finite known-conflict margins and role-gain deltas are descriptive TRAIN window-chain '
            'diagnostics, not full-grammar regret, biological false positives, gradient-matched causality, '
            'raw-FASTA accuracy, independent generalization or novelty.', **D.SCOPE))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--checks-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.checks_job)
