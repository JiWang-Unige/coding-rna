#!/usr/bin/env python3
"""A35 descriptive components of the closed A34 potential, never a new score."""
from collections import Counter, defaultdict
import json
import math
import time
import torch
import diagnose_missing_path_emission_train as D

OUT = D.ROOT/'outputs/M28-A35-SAVED-ROLE-CONTRIBUTIONS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a35_saved_role_contributions_protocol.md'
GROUPS = dict(coding=(1, 2, 3), intron=(4, 5, 6), donor=(7, 8, 9),
              acceptor=(10, 11, 12), start=(13,), stop=(14,))
STRATA = D.STATES+('natural_hit_control',)


@torch.inference_mode()
def components(logits, chain):
    roles = D.R.candidate_roles(chain, len(logits))
    positions = torch.arange(len(logits))
    values = logits.double()[positions, roles]-logits[:, 0].double()
    role_sums = torch.bincount(roles, weights=values, minlength=15).tolist()
    role_bases = torch.bincount(roles, minlength=15).tolist()
    groups = {}
    for name, indices in GROUPS.items():
        bases = sum(role_bases[i] for i in indices); gain = sum(role_sums[i] for i in indices)
        groups[name] = dict(role_bases=bases, gain=gain, gain_per_role_base=gain/bases if bases else None)
    negative = [name for name in GROUPS if groups[name]['role_bases'] and groups[name]['gain'] < 0]
    dominant = min(negative, key=lambda name: groups[name]['gain']) if negative else None
    return dict(groups=groups, fine_role_gain=role_sums[1:], fine_role_bases=role_bases[1:],
        component_sum=sum(g['gain'] for g in groups.values()), dominant_negative_group=dominant)


def statistics(rows):
    groups = {}
    for name in GROUPS:
        bases = sum(r['components']['groups'][name]['role_bases'] for r in rows)
        gain = sum(r['components']['groups'][name]['gain'] for r in rows)
        groups[name] = dict(role_bases=bases, gain_sum=gain,
            pooled_gain_per_role_base=gain/bases if bases else None,
            pairs_present=sum(r['components']['groups'][name]['role_bases'] > 0 for r in rows),
            pairs_negative=sum(r['components']['groups'][name]['gain'] < 0 for r in rows))
    dominant = Counter(r['components']['dominant_negative_group'] for r in rows)
    return dict(reference_window_pairs=len(rows), groups=groups,
        dominant_negative_group_counts={name: dominant[name] for name in GROUPS},
        no_negative_group=dominant[None])


def aggregate(rows):
    return dict(all=statistics(rows), by_original_A34_state={s: statistics([r for r in rows
        if r['original_A34_state'] == s]) for s in STRATA},
        per_species={species: {s: statistics([r for r in rows
            if r['species'] == species and r['original_A34_state'] == s]) for s in STRATA}
            for species in D.R.C.R.ALLOW},
        nonpositive_by_exon_class={c: statistics([r for r in rows if r['original_A34_state'] == D.STATES[0]
            and r['exon_class'] == c]) for c in ('one', 'two_to_four', 'five_or_more')})


@torch.inference_mode()
def main():
    began = time.perf_counter(); D.require_jobs(['13403916', '13403920'])
    source = D.R.completed(D.OUT, 13403916); checked = D.R.completed(D.CHECK, 13403920)
    if (source['protocol'] != D.PROTOCOL or checked['protocol'] != D.PROTOCOL or
            source['phase'] != 'score' or checked['phase'] != 'check' or
            source['counts'] != checked['counts'] or source['windows'] != 396 or
            source['reference_window_pairs'] != 860 or not source['model_weights_unchanged'] or
            not checked['actual_saved_reference_identity_and_direct_scores_verified'] or
            not checked['actual_single_and_joint_reference_assisted_selection_recomputed'] or
            not checked['actual_saved_aggregate_and_original_denominators_verified']):
        raise ValueError('Actual closed full original A34 scoring and independent check required')
    refs = D.R.C.R.rows(D.OUT/'references.jsonl'); windows = D.R.C.R.rows(D.OUT/'windows.jsonl')
    grouped = defaultdict(list)
    for ref in refs:
        grouped[ref['window_id']].append(ref)
    if (len(refs) != 860 or len(windows) != 396 or len(grouped) != 396 or
            [w['window']['id'] for w in windows] != list(grouped) or
            len({(r['window_id'], r['original_candidate_index']) for r in refs}) != 860):
        raise ValueError('Original full A34 reference/window identity population required')
    rows = []; max_error = 0.
    with (OUT/'references.jsonl').open('x') as handle:
        for ordinal, win in enumerate(windows):
            w = win['window']; ident = w['id']
            expected = D.OUT/'logits'/('%04d.pt' % ordinal)
            if (w['split'] != 'train' or w['seqid'] not in D.R.C.R.ALLOW[w['species']] or
                    win['logits_artifact'] != str(expected.relative_to(D.ROOT))):
                raise ValueError('Original TRAIN-only saved logit path required')
            state = torch.load(expected, map_location='cpu', weights_only=True); logits = state['logits']
            if state['window'] != w or logits.dtype != torch.float32 or tuple(logits.shape) != (24576, 15):
                raise ValueError('Unchanged valid FP32 saved emission geometry required')
            for ref in grouped[ident]:
                parts = components(logits, ref['chain']); old = ref['reference_gain']
                if not math.isclose(parts['component_sum'], old, rel_tol=1e-12, abs_tol=1e-8):
                    raise ValueError('Original A34 gain not reproduced by full role decomposition')
                max_error = max(max_error, abs(parts['component_sum']-old))
                row = {k: ref[k] for k in ('window_id', 'original_candidate_index', 'species', 'chain',
                    'exon_class', 'cohort', 'original_local_mechanism')}
                row.update(original_A34_state=ref.get('reference_assisted_single_state', 'natural_hit_control'),
                           original_A34_gain=old, components=parts)
                rows.append(row); handle.write(json.dumps(row, allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0 or ordinal+1 == len(windows):
                handle.flush(); print(json.dumps(dict(windows=ordinal+1, pairs=len(rows),
                    seconds=time.perf_counter()-began)), flush=True)
    expected_strata = {**source['counts']['missing']['states'], 'natural_hit_control': 324}
    actual_strata = Counter(r['original_A34_state'] for r in rows)
    if (len(rows) != 860 or {s: actual_strata[s] for s in STRATA} != expected_strata or
            D.R.C.R.rows(OUT/'references.jsonl') != rows):
        raise ValueError('Actual saved complete records or original A34 strata differ')
    result = aggregate(rows)
    D.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A35-SAVED-ROLE-CONTRIBUTIONS',
        protocol=PROTOCOL, source_jobs=dict(score='13403916', CPU='13403920'), windows=396,
        reference_window_pairs=860, role_groups=GROUPS, original_A34_counts=source['counts'],
        contribution_statistics=result, max_component_sum_vs_original_error=max_error,
        actual_original_scores_and_strata_reproduced=True, actual_saved_complete_records_verified=True,
        original_model_scores_weights_threshold_and_selection_unchanged=True,
        model_or_backbone_forward=False, feature_cache_or_catalog_opened=False,
        optimizer_steps=0, new_parameters=0, new_candidate_search=False, new_R4_labelling=False,
        new_caller_or_rescored_predictions=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, seconds=time.perf_counter()-began,
        limitation='Descriptive original in-sample role contributions, not a modified score, '
            'causal role ablation, biological mechanism, free prediction, performance gate or generalization.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
