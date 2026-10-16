#!/usr/bin/env python3
"""A34 old TRAIN reference-assisted scoring, never free candidate generation."""
import argparse
from collections import Counter, defaultdict
import json
import math
import subprocess
import time
import torch
import run_local_path_caller as R

ROOT = R.ROOT
OUT = ROOT/'outputs/M28-A34-MISSING-PATH-EMISSION-DIAGNOSTIC'
CHECK = ROOT/'outputs/M28-A34-MISSING-PATH-EMISSION-CHECK'
A27 = ROOT/'outputs/M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC'
A24 = ROOT/'outputs/M28-A24-COMPLETED-PREFIX-COVERAGE'
PROTOCOL = 'reports/M28-METHOD-RESTART/a34_missing_path_emission_protocol.md'
STATES = ('nonpositive_gain', 'positive_not_selected',
          'positive_selected_without_original_correct_displacement',
          'positive_selected_with_original_correct_displacement')


def key(chain):
    return tuple(map(tuple, chain))


def require_jobs(jobs):
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {f[0].strip(): [x.strip() for x in f[1:3]]
              for line in result.stdout.splitlines() if len(f := line.split('|')) >= 3}
    if any(actual.get(j) != ['COMPLETED', '0:0'] for j in jobs):
        raise ValueError('Closed original source jobs required: '+str(actual))


def sources():
    require_jobs(['13391716', '13392411', '13398521', '13398524'])
    a27 = R.completed(A27, 13392411); a24 = R.completed(A24, 13391716)
    prior = R.completed(R.READOUT, 13398521); checked = R.completed(R.CHECK, 13398524)
    if (a27['generator_checkpoint'] != R.C.GENERATOR or a27['windows'] != 396 or
            a27['missing_reference_window_pairs'] != 536 or a27['control_reference_window_pairs'] != 324 or
            not a27['actual_original_missing_and_control_identity_verified'] or
            not a27['actual_saved_profiles_geometry_components_first_exclusions_and_summaries_reconciled'] or
            not a27['all_saved_natural_beam_controls_pass_necessary_local_conditions'] or
            a24['pooled']['reference_window_pairs'] != 2916 or
            a24['pooled']['missing_pairs'] != 536 or
            prior['configuration'] != R.CONFIG or checked['configuration'] != R.CONFIG or
            prior['windows'] != 1491 or checked['windows'] != 1491 or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] or
            not checked['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'] or
            prior['selected_counts_per_species'] != checked['selected_counts_per_species'] or
            prior['eligible_reference_pairs'] != checked['eligible_reference_pairs']):
        raise ValueError('Matching full original A27 and A31 evidence required')
    ids = list(dict.fromkeys(d['window_id'] for d in R.C.R.rows(R.Q.DRAWS)))
    natural_rows = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    label_rows = R.C.R.rows(R.READOUT/'windows.jsonl')
    free_rows = R.C.R.rows(R.C.FIT/'free_pool.jsonl')
    if (len(ids) != 1491 or any([x['window']['id'] for x in rows] != ids
                              for rows in (natural_rows, label_rows, free_rows)) or
            sum(len(x['chains']) for x in free_rows) != 16321):
        raise ValueError('Original complete ordered TRAIN population required')
    natural = {x['window']['id']: x for x in natural_rows}
    labels = {x['window']['id']: x for x in label_rows}
    free = {x['window']['id']: x for x in free_rows}
    original = defaultdict(list)
    for ref in R.C.R.rows(A24/'references.jsonl'):
        original[ref['window_id']].append(ref)
    refs = R.C.R.rows(A27/'references.jsonl')
    expected = []
    for ident in ids:
        missing = [r for r in original[ident] if r['category'] != 'natural_hit']
        if not missing:
            continue
        geometry = {key(c): i for i, c in enumerate(free[ident]['chains'])}
        controls = [r for r in original[ident] if r['category'] == 'natural_hit' and
                    'beam4' in free[ident]['decoder_sources'][geometry[key(r['chain'])]]]
        expected.extend(missing+controls[:1])
    identity = ('window_id', 'original_candidate_index', 'species', 'chain', 'origin', 'exon_class', 'category')
    if (len(refs) != 860 or len(expected) != 860 or
            any(any(a[k] != b[k] for k in identity) for a, b in zip(refs, expected)) or
            len({(r['window_id'], key(r['chain'])) for r in refs}) != 860):
        raise ValueError('All original missing pairs and metadata-first controls required')
    grouped = defaultdict(list)
    for ref in refs:
        ident = ref['window_id']; w = natural[ident]['window']; pool = natural[ident]['variants']['A31_A3']
        if (w['split'] != 'train' or w['seqid'] not in R.C.R.ALLOW[w['species']] or
                w['valid_bases'] != 24576 or w['width'] != 24576 or
                any(w[k] != labels[ident]['window'][k] or w[k] != free[ident]['window'][k]
                    for k in R.C.R.GEOMETRY) or natural[ident]['reference_used_in_scoring'] or
                free[ident]['reference_used_in_pool_generation'] or pool['chains'] != free[ident]['chains'] or
                not ref['grammar_supported']):
            raise ValueError('Unchanged original TRAIN-only pool and supported oriented reference required')
        position = {key(c): i for i, c in enumerate(pool['chains'])}.get(key(ref['chain']))
        if ref['cohort'] == 'missing':
            if position is not None:
                raise ValueError('Original missing geometry was already in natural pool')
        elif ref['cohort'] != 'natural_hit_control' or position is None or labels[ident]['labels']['A31_A3'][position] != 1:
            raise ValueError('Original natural control must be a saved correct pool member')
        grouped[ident].append(ref)
    missing = [r for r in refs if r['cohort'] == 'missing']
    if (len(grouped) != 396 or len(missing) != 536 or
            Counter(r['mechanism'] for r in missing) != Counter({
                'all_queries_locally_excluded': 315, 'some_query_passes_local_necessary_conditions': 221}) or
            Counter(r['species'] for r in missing) != Counter(arabidopsis_thaliana=421, oryza_sativa=115)):
        raise ValueError('Original diagnostic denominators changed')
    index = {r['window_id']: r for r in R.Q.train_records(R.Q.FEATURES/'index.jsonl')
             if r['window_id'] in grouped}
    if set(index) != set(grouped):
        raise ValueError('All original TRAIN-only cached feature identities required')
    return grouped, natural, labels, index


def assisted_selection(pool, labels, additions):
    """Fixed scores/DP in an explicitly reference-assisted diagnostic pool."""
    if (len(pool['chains']) != len(pool['scores']) or len(labels) != len(pool['chains']) or
            R.select_nonoverlapping(pool['chains'], pool['scores']) != pool['selected']):
        raise ValueError('Saved original selection/label identity differs')
    entries = {key(c): (s, y, False) for c, s, y in zip(pool['chains'], pool['scores'], labels)}
    if len(entries) != len(pool['chains']) or list(entries) != sorted(entries):
        raise ValueError('Original unique canonical geometry order required')
    for chain, gain in additions:
        k = key(chain)
        if k in entries or not math.isfinite(gain):
            raise ValueError('Only distinct missing finite-scored reference chains may be added')
        entries[k] = (gain, 1, True)
    chains = sorted(entries)
    chosen = R.select_nonoverlapping(chains, [entries[k][0] for k in chains])
    selected = [chains[i] for i in chosen]
    old_selected = {key(pool['chains'][i]) for i in pool['selected']}
    dropped = old_selected-set(selected)
    removed = {str(y): sum(entries[k][1] == y for k in dropped) for y in (1, 0, -1)}
    counts = {str(y): sum(entries[k][1] == y for k in selected) for y in (1, 0, -1)}
    return dict(selected_chains=[[list(e) for e in c] for c in selected], selected_label_counts=counts,
        injected_references_selected=sum(entries[k][2] for k in selected),
        original_selected_dropped_label_counts=removed)


def evaluate_window(refs, gains, pool, labels):
    if len(refs) != len(gains):
        raise ValueError('One original score per ordered reference required')
    results = []; missing = []
    geometry = {key(c): i for i, c in enumerate(pool['chains'])}
    for ref, gain in zip(refs, gains):
        if ref['cohort'] == 'missing':
            single = assisted_selection(pool, labels, [(ref['chain'], gain)])
            picked = key(ref['chain']) in {key(c) for c in single['selected_chains']}
            if gain <= 0:
                if picked:
                    raise ValueError('Original null0 cannot select a nonpositive reference')
                state = STATES[0]
            elif not picked:
                state = STATES[1]
            else:
                state = STATES[3 if single['original_selected_dropped_label_counts']['1'] else 2]
            result = dict(reference_gain=gain, reference_assisted_single_state=state,
                          single_reference_assisted_selection=single)
            missing.append((ref['chain'], gain))
        else:
            i = geometry[key(ref['chain'])]
            if not math.isclose(gain, pool['scores'][i], rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError('Natural control direct score differs from original saved pool')
            result = dict(reference_gain=gain, original_control_selected=i in pool['selected'])
        results.append(dict(window_id=ref['window_id'], original_candidate_index=ref['original_candidate_index'],
            species=ref['species'], chain=ref['chain'], exon_class=ref['exon_class'],
            cohort=ref['cohort'], original_local_mechanism=ref['mechanism'], **result))
    joint = assisted_selection(pool, labels, missing)
    baseline = {str(y): sum(labels[i] == y for i in pool['selected']) for y in (1, 0, -1)}
    joint['original_selected_label_counts'] = baseline
    joint['selected_label_count_changes'] = {str(y): joint['selected_label_counts'][str(y)]-baseline[str(y)]
                                            for y in (1, 0, -1)}
    joint['injected_missing_reference_geometries'] = len(missing)
    return results, joint


def aggregate(records, windows):
    missing = [r for r in records if r['cohort'] == 'missing']
    def counts(rows):
        c = Counter(r['reference_assisted_single_state'] for r in rows)
        return dict(reference_window_pairs=len(rows), states={s: c[s] for s in STATES})
    controls = [r for r in records if r['cohort'] == 'natural_hit_control']
    return dict(missing=counts(missing), per_species={s: counts([r for r in missing if r['species'] == s])
        for s in R.C.R.ALLOW}, by_original_local_mechanism={m: counts([r for r in missing
        if r['original_local_mechanism'] == m]) for m in (
            'all_queries_locally_excluded', 'some_query_passes_local_necessary_conditions')},
        by_exon_class={c: counts([r for r in missing if r['exon_class'] == c])
            for c in ('one', 'two_to_four', 'five_or_more')},
        natural_hit_controls=dict(reference_window_pairs=len(controls),
            selected=sum(r['original_control_selected'] for r in controls),
            positive_gain=sum(r['reference_gain'] > 0 for r in controls)),
        joint_reference_assisted=dict(windows=len(windows),
            injected_missing_reference_geometries=sum(w['joint']['injected_missing_reference_geometries'] for w in windows),
            injected_references_selected=sum(w['joint']['injected_references_selected'] for w in windows),
            original_selected_dropped_label_counts={str(y): sum(w['joint']['original_selected_dropped_label_counts'][str(y)]
                for w in windows) for y in (1, 0, -1)},
            selected_label_count_changes={str(y): sum(w['joint']['selected_label_count_changes'][str(y)]
                for w in windows) for y in (1, 0, -1)}))


def write(out, **values):
    R.C.R.write_summary(out, dict(completed=True, experiment='M28-A34-MISSING-PATH-EMISSION-DIAGNOSTIC',
        protocol=PROTOCOL, generator_checkpoint=R.C.GENERATOR,
        source_jobs=dict(A24='13391716', A27='13392411', A31='13398521', A31_CPU='13398524'),
        windows=396, reference_window_pairs=860, missing_pairs=536, control_pairs=324,
        score=R.CONFIG['score'], null=0., all_role_weights=1.,
        reference_assisted=True, reference_injected_diagnostic_geometries=536,
        reference_injected_actual_caller_candidates=0, original_natural_predictions_unchanged=True,
        optimizer_steps=0, new_parameters=0, new_candidate_search=False, new_backbone_forward=False,
        new_R4_labelling=False, catalog_opened=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, A33_search_checker_or_stop_condition_changed=False,
        limitation='Old in-sample TRAIN reference/window diagnostic, not free predictions, '
            'gene F1, jointly summable individual recoveries, independent generalization or a model PASS.', **values))


@torch.inference_mode()
def score():
    began = time.perf_counter(); grouped, natural, labels, index = sources()
    model = R.C.load_generator(); before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649457:
        raise ValueError('Original frozen A3 capacity required')
    (OUT/'logits').mkdir(); records = []; windows = []; max_error = 0.
    torch.cuda.reset_peak_memory_stats()
    with (OUT/'references.jsonl').open('x') as ref_file, (OUT/'windows.jsonl').open('x') as win_file:
        for ordinal, (ident, refs) in enumerate(grouped.items()):
            w = natural[ident]['window']; pool = natural[ident]['variants']['A31_A3']
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            R.Q.validate_cached(cached, w)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]; logits = model.segmentation(hidden)
            if logits.dtype != torch.float32:
                raise ValueError('Original FP32 forward required')
            old = R.local_path_scores(logits, pool['chains'], w['valid_bases']).cpu()
            expected = torch.tensor(pool['scores'], dtype=torch.float64)
            torch.testing.assert_close(old, expected, rtol=1e-12, atol=1e-8)
            max_error = max(max_error, max((abs(float(a-b)) for a, b in zip(old, expected)), default=0.))
            gains = R.local_path_scores(logits, [r['chain'] for r in refs], w['valid_bases']).cpu().tolist()
            result, joint = evaluate_window(refs, gains, pool, labels[ident]['labels']['A31_A3'])
            path = OUT/'logits'/('%04d.pt' % ordinal)
            torch.save(dict(window=w, logits=logits[:w['valid_bases']].detach().cpu().clone()), path)
            window = dict(window=w, logits_artifact=str(path.relative_to(ROOT)), joint=joint)
            for row in result:
                ref_file.write(json.dumps(row, allow_nan=False)+'\n')
            win_file.write(json.dumps(window, allow_nan=False)+'\n')
            records.extend(result); windows.append(window)
            if (ordinal+1) % 32 == 0 or ordinal+1 == len(grouped):
                ref_file.flush(); win_file.flush()
                print(json.dumps(dict(windows=ordinal+1, reference_pairs=len(records),
                                      seconds=time.perf_counter()-began)), flush=True)
            del cached, features, dna, hidden, logits
    if (len(records) != 860 or len(windows) != 396 or
            any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Incomplete original coverage or frozen generator changed')
    write(OUT, phase='score', counts=aggregate(records, windows), model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True, original_A31_A3_scores_reproduced=True,
        max_original_score_error=max_error, model_forward=True, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated())


def check(source_job):
    began = time.perf_counter(); require_jobs([source_job]); report = R.completed(OUT, source_job)
    if report['phase'] != 'score' or report['protocol'] != PROTOCOL or not report['model_weights_unchanged']:
        raise ValueError('Matching completed actual A34 score source required')
    grouped, natural, labels, _ = sources()
    saved = R.C.R.rows(OUT/'references.jsonl'); windows = R.C.R.rows(OUT/'windows.jsonl')
    if (len(saved) != 860 or len(windows) != 396 or
            [w['window']['id'] for w in windows] != list(grouped)):
        raise ValueError('Actual complete original saved output coverage required')
    actual = []; checked_windows = []; max_error = 0.
    for ordinal, win in enumerate(windows):
        ident = win['window']['id']; w = natural[ident]['window']; refs = grouped[ident]
        if win['window'] != w or win['logits_artifact'] != str((OUT/'logits'/('%04d.pt' % ordinal)).relative_to(ROOT)):
            raise ValueError('Actual saved logit geometry/identity differs')
        state = torch.load(ROOT/win['logits_artifact'], map_location='cpu', weights_only=True)
        logits = state['logits']
        if state['window'] != w or logits.dtype != torch.float32 or tuple(logits.shape) != (24576, 15):
            raise ValueError('Original valid oriented FP32 saved emissions required')
        logp = logits.double().log_softmax(-1); positions = torch.arange(len(logp))
        pool = natural[ident]['variants']['A31_A3']; values = []
        for chain in pool['chains']+[r['chain'] for r in refs]:
            roles = R.candidate_roles(chain, len(logp))
            values.append(float((logp[positions, roles]-logp[:, 0]).sum()))
        old = values[:len(pool['chains'])]; gains = values[len(pool['chains']):]
        torch.testing.assert_close(torch.tensor(old, dtype=torch.float64),
            torch.tensor(pool['scores'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
        result, joint = evaluate_window(refs, gains, pool, labels[ident]['labels']['A31_A3'])
        expected = saved[len(actual):len(actual)+len(result)]
        for row, prior in zip(result, expected):
            error = abs(row['reference_gain']-prior['reference_gain']); max_error = max(max_error, error)
            if not math.isclose(row['reference_gain'], prior['reference_gain'], rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError('Direct saved log-probability difference differs from scored reference')
            # Classification/DP must reproduce the actual saved FP64 potential,
            # not flip a boundary on a different algebraically equivalent sum.
        replay, replay_joint = evaluate_window(refs, [r['reference_gain'] for r in expected],
            pool, labels[ident]['labels']['A31_A3'])
        if replay != expected or replay_joint != win['joint']:
            raise ValueError('Actual saved scores, reference classes or assisted interval selections differ')
        # The independent direct sum must also yield the same selection/category.
        for row, prior in zip(result, expected):
            row['reference_gain'] = prior['reference_gain']
        if result != expected or joint != win['joint']:
            raise ValueError('Independent direct scores changed a null/selection/category decision')
        actual.extend(replay); checked_windows.append(win)
        if (ordinal+1) % 64 == 0 or ordinal+1 == len(windows):
            print(json.dumps(dict(CPU_windows=ordinal+1, seconds=time.perf_counter()-began)), flush=True)
    counts = aggregate(actual, checked_windows)
    if counts != report['counts']:
        raise ValueError('Actual saved independent aggregate differs')
    write(CHECK, phase='check', score_job=source_job, counts=counts,
        actual_saved_reference_identity_and_direct_scores_verified=True,
        actual_single_and_joint_reference_assisted_selection_recomputed=True,
        actual_saved_aggregate_and_original_denominators_verified=True,
        max_direct_reference_score_error=max_error, model_forward=False,
        seconds=time.perf_counter()-began)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('score', 'check'), required=True)
    parser.add_argument('--source-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'score':
        score()
    elif args.source_job:
        check(args.source_job)
    else:
        parser.error('Actual --source-job required for CPU check')
