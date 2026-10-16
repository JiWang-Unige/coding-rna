#!/usr/bin/env python3
"""Read-only A23 actual exposure, saved chains, TRAIN labels and fixed gate."""
import argparse
from collections import Counter
import json
import math
import torch
import run_query_relative_boundary as R
from src.m28.query_model import unmatched_query_policy

OUT = R.ROOT/'outputs/M28-A23-RELATIVE-BOUNDARY-RESULT-CHECK'


def main(fit_job, readout_job):
    for folder, job in ((R.FIT, fit_job), (R.READOUT, readout_job)):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    fitted = json.loads((R.FIT/'summary.json').read_text()); R.validate_fit(fitted)
    report = json.loads((R.READOUT/'summary.json').read_text())
    assert report['completed'] and report['arm'] == 'A23' and report['configuration'] == R.CONFIG and report['protocol'] == R.PROTOCOL
    assert report['checkpoint'] == fitted['final_checkpoint'] and report['optimizer_steps_this_run'] == 0
    assert report['new_candidate_search'] and report['original_A3_actual_natural_pool_labels_recomputed']
    assert not report['reference_used_in_pool_generation'] and report['reference_injected_candidates'] == 0
    assert not report['new_backbone_forward'] and not report['DEV_forward_used'] and not report['test_setaria_used']
    state = torch.load(R.ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R.validate_fit(state, state=True)
    assert all(bool(torch.isfinite(v).all()) for v in state['model'].values())
    model = R.new_model('cpu'); model.load_state_dict(state['model'], strict=True)
    assert sum(p.numel() for p in model.parameters()) == fitted['trainable_parameters'] == 649617
    fresh = R.new_model('cpu')
    changes = {k: float((v-fresh.query.relative_boundary.state_dict()[k]).double().norm())
               for k, v in model.query.relative_boundary.state_dict().items()}
    for k, value in changes.items():
        assert math.isclose(value, fitted['relative_parameter_change_norms'][k], rel_tol=1e-12, abs_tol=1e-12)
    assert changes['0.weight'] > 0 and changes['2.weight'] > 0
    del model, fresh, state
    for step in range(256, 4609, 256):
        assert (R.FIT/('step_'+str(step).zfill(6)+'.pt')).stat().st_size > 0
    draws, index, windows, catalogs = R.Q.load_train_inputs()
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    assert len(draws) == 1536
    free = R.C.R.rows(R.READOUT/'free_pool.jsonl'); records = R.C.R.rows(R.READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in free] == [r['window']['id'] for r in records] == ids
    assert len(free) == len(records) == report['windows'] == 1491
    baseline_free = {r['window']['id']: r for r in R.C.R.rows(R.C.FIT/'free_pool.jsonl')}
    old = {r['window']['id']: r for r in R.C.R.rows(R.D.OUT/'windows.jsonl')}
    actual = {s: {a: Counter() for a in ('A3', 'A23')} for s in R.A.ALLOW}
    supervision = {}; candidate_pairs = 0
    for raw, rec in zip(free, records):
        ident = raw['window']['id']; w = windows[ident]
        assert rec['window'] == raw['window'] and not raw['reference_used_in_pool_generation']
        assert w['split'] == 'train' and w['seqid'] in R.C.R.ALLOW[w['species']]
        assert all(raw['window'][k] == w[k] == old[ident]['window'][k] for k in R.C.R.GEOMETRY)
        assert [r['chain'] for r in rec['candidates']] == raw['chains']
        assert len({tuple(map(tuple, chain)) for chain in raw['chains']}) == len(raw['chains'])
        assert len(raw['decoder_sources']) == len(raw['chains'])
        cached = torch.load(R.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        R.Q.validate_cached(cached, w)
        targets = R.build_targets(w, catalogs[w['species']], cached['sequence'])
        assert R.candidate_targets(raw['chains'], targets) == [r['label'] for r in rec['candidates']]
        baseline = [r for r in old[ident]['candidates'] if r['origin'] == 'free']
        assert [r['chain'] for r in baseline] == baseline_free[ident]['chains']
        assert R.candidate_targets([r['chain'] for r in baseline], targets) == [r['label'] for r in baseline]
        values = dict(A3=R.counts(baseline), A23=R.counts(rec['candidates'])); assert values == rec['natural_counts']
        for arm in ('A3', 'A23'):
            actual[w['species']][arm].update(values[arm])
        n = len(targets['positive_chains']); policy = unmatched_query_policy(targets, w['valid_bases'])
        supervision[ident] = dict(counts={'1': n, '0': 16-n if policy['unmatched_label'] == 0 else 0,
                                          '-1': 16-n if policy['unmatched_label'] == -1 else 0},
                                 protection_reasons=list(policy['protection_reasons']))
        candidate_pairs += len(raw['chains']); del cached, targets
    actual = {s: {a: dict(c) for a, c in arms.items()} for s, arms in actual.items()}
    assert actual == report['natural_counts_per_species']
    checks = R.gate(actual)
    assert checks == report['predeclared_TRAIN_generator_gate'] and all(checks.values()) == report['TRAIN_generator_supported']
    training = R.C.R.rows(R.FIT/'training.jsonl'); assert len(training) == 4608
    counts = [{s: Counter() for s in R.A.ALLOW} for _ in range(3)]; losses = [Counter() for _ in range(3)]
    gradient_steps = Counter()
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; ident = draw['window_id']; species = windows[ident]['species']
        assert (row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) == (i+1, epoch+1, ordinal, ident, species)
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        assert row['query_label_counts'] == supervision[ident]['counts'] and row['protection_reasons'] == supervision[ident]['protection_reasons']
        assignment = row['gene_to_query']; assert len(assignment) == len(set(assignment)) == row['query_label_counts']['1']
        assert all(0 <= q < 16 for q in assignment)
        relative = row['relative_layer_gradient_norms_before_clip']; assert set(relative) == {'input', 'output'}
        assert all(math.isfinite(v) for v in (*row['losses'].values(), row['weighted_total'], row['gradient_norm_before_clip'], *relative.values()))
        assert all(v >= 0 for v in relative.values()); gradient_steps.update({k: int(v > 0) for k, v in relative.items()})
        expected = row['importance_weight']*sum(row['losses'].values())
        assert abs(expected-row['weighted_total']) <= 1e-5*(1+abs(row['weighted_total']))
        counts[epoch][species].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
        counts[epoch][species].update(row['protection_reasons']); losses[epoch].update(dict(row['losses'], weighted_total=row['weighted_total']))
    assert dict(gradient_steps) == fitted['relative_layer_nonzero_gradient_steps'] and all(v > 0 for v in gradient_steps.values())
    original = json.loads((R.A.STAGE_DIRS['fit']/'summary.json').read_text()); R.A.check_last_identity(original)
    assert len(fitted['epochs']) == len(original['epochs']) == 3
    for i in range(3):
        actual_supervision = {s: dict(c) for s, c in counts[i].items()}
        assert actual_supervision == fitted['epochs'][i]['supervision_by_species'] == original['epochs'][i]['supervision_by_species']
        assert fitted['epochs'][i]['steps'] == 1536
        assert {k: v/1536 for k, v in losses[i].items()} == fitted['epochs'][i]['mean_losses']
    result = dict(completed=True, arm='A23', source_jobs=dict(fit=fit_job, readout=readout_job), windows=1491,
        actual_4608_draw_weight_and_A3_supervision_identity_verified=True,
        actual_new_and_original_natural_chain_labels_recomputed=True,
        actual_saved_pool_identity_and_three_fixed_gate_checks_verified=True,
        last_model_structure_and_all_tensors_finite=True, trainable_parameters=649617,
        actual_relative_layers_gradients_and_weight_changes_verified=True,
        reference_injected_candidates=0, candidate_window_pairs=candidate_pairs,
        natural_counts_per_species=actual, original_frozen_TRAIN_generator_gate=checks,
        TRAIN_generator_supported=all(checks.values()), optimizer_steps_this_run=0,
        new_head_or_backbone_forward=False, new_candidate_search=False, DEV_forward_used=False, test_setaria_used=False,
        limitation='In-sample natural candidate coverage only; no selected gene F1 or independent generalization.')
    with (OUT/'summary.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-job', required=True); parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.fit_job, args.readout_job)
