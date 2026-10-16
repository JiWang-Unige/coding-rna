#!/usr/bin/env python3
"""CPU-only A30 exposure, safe-negative records and five frozen TRAIN gates."""
import argparse
from collections import Counter
import json
import math
import torch
import run_query_search_retention as R
from src.m28.query_model import unmatched_query_policy
from query_search_retention_training import loss_nodes, classify, compatible

OUT = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-RESULT-CHECK'


def main(fit_job, readout_job):
    for folder, job in ((R.FIT, fit_job), (R.READOUT, readout_job)):
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    fitted = json.loads((R.FIT/'summary.json').read_text()); R.validate_fit(fitted)
    report = json.loads((R.READOUT/'summary.json').read_text())
    assert report['completed'] and report['arm'] == 'A30' and report['configuration'] == R.CONFIG and report['protocol'] == R.PROTOCOL
    assert report['checkpoint'] == fitted['final_checkpoint'] and report['optimizer_steps_this_run'] == 0
    assert report['new_candidate_search'] and report['original_A3_actual_natural_pool_labels_recomputed']
    assert not report['reference_used_in_pool_generation'] and report['reference_injected_candidates'] == 0
    assert not report['new_backbone_forward'] and not report['DEV_forward_used'] and not report['test_setaria_used']
    state = torch.load(R.ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R.validate_fit(state, state=True)
    assert all(bool(torch.isfinite(v).all()) for v in state['model'].values())
    model = R.new_model('cpu'); model.load_state_dict(state['model'], strict=True)
    assert sum(p.numel() for p in model.parameters()) == fitted['trainable_parameters'] == 649457
    del model, state
    for step in range(256, 4609, 256):
        assert (R.FIT/('step_'+str(step).zfill(6)+'.pt')).stat().st_size > 0
    draws, index, windows, catalogs = R.Q.load_train_inputs()
    ids = list(dict.fromkeys(d['window_id'] for d in draws)); assert len(draws) == 1536
    free = R.C.R.rows(R.READOUT/'free_pool.jsonl'); records = R.C.R.rows(R.READOUT/'windows.jsonl')
    assert [r['window']['id'] for r in free] == [r['window']['id'] for r in records] == ids
    assert len(free) == len(records) == report['windows'] == 1491
    baseline_free = {r['window']['id']: r for r in R.C.R.rows(R.C.FIT/'free_pool.jsonl')}
    old = {r['window']['id']: r for r in R.C.R.rows(R.D.OUT/'windows.jsonl')}
    actual = {s: {a: Counter() for a in ('A3', 'A30')} for s in R.A.ALLOW}
    supervision = {}; window_targets = {}; candidate_pairs = 0
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
        values = dict(A3=R.counts(baseline), A30=R.counts(rec['candidates'])); assert values == rec['natural_counts']
        for arm in ('A3', 'A30'):
            actual[w['species']][arm].update(values[arm])
        n = len(targets['positive_chains']); policy = unmatched_query_policy(targets, w['valid_bases'])
        supervision[ident] = dict(counts={'1': n, '0': 16-n if policy['unmatched_label'] == 0 else 0,
                                          '-1': 16-n if policy['unmatched_label'] == -1 else 0},
                                 protection_reasons=list(policy['protection_reasons']))
        window_targets[ident] = {k: targets[k] for k in
            ('positive_chains', 'protected_chains', 'unknown_chain_spans', 'callable_bases')}
        candidate_pairs += len(raw['chains']); del cached, targets
    actual = {s: {a: dict(c) for a, c in arms.items()} for s, arms in actual.items()}
    assert actual == report['natural_counts_per_species']
    # Confirm the declared baseline counts before evaluating the new gates.
    expected_baselines = {'arabidopsis_thaliana': (1751, 6305, 604), 'oryza_sativa': (629, 6261, 184)}
    for species, expected in expected_baselines.items():
        assert tuple(actual[species]['A3'][k] for k in ('correct', 'known_nonmatch', 'five_or_more_correct')) == expected
    checks = R.gate(actual)
    assert len(checks) == 5 and checks == report['predeclared_TRAIN_generator_gate']
    assert all(checks.values()) == report['TRAIN_generator_supported']
    training = R.C.R.rows(R.FIT/'training.jsonl'); assert len(training) == 4608
    counts = [{s: Counter() for s in R.A.ALLOW} for _ in range(3)]; losses = [Counter() for _ in range(3)]
    margin_steps = Counter(); dispositions = Counter(); qualified_pairs = 0
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; ident = draw['window_id']; species = windows[ident]['species']
        assert (row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) == (i+1, epoch+1, ordinal, ident, species)
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        assert row['query_label_counts'] == supervision[ident]['counts'] and row['protection_reasons'] == supervision[ident]['protection_reasons']
        assignment = row['gene_to_query']; assert len(assignment) == len(set(assignment)) == row['query_label_counts']['1']
        assert all(0 <= q < 16 for q in assignment)
        targets = window_targets[ident]
        assigned = {q: tuple(map(tuple, c)) for q, c in zip(assignment, targets['positive_chains'])}
        collection = row['collection_records']
        assert row['positive_queries'] == len(collection) == len(assigned)
        assert [r['query_index'] for r in collection] == sorted(assigned) and row['reference_injected_candidates'] == 0
        terms = []
        for record in collection:
            reference = tuple(map(tuple, record['chain']))
            assert reference == assigned[record['query_index']]
            disposition = record['supervision_disposition']; dispositions[disposition] += 1
            if not record['raw_grammar_supported']:
                assert disposition == 'original_teacher_only' and record['unsupported_reason']
                terms.append(0.); continue
            if record['returned']:
                assert disposition == 'already_returned' and record['first_loss'] is None
                terms.append(0.); continue
            target, cutoff = loss_nodes(record['first_loss'])
            assert json.loads(json.dumps(target)) == record['target_node']
            assert json.loads(json.dumps(cutoff)) == record['cutoff_node']
            assert compatible(target, reference)
            annotation = classify(cutoff, targets, reference)
            assert annotation == record['cutoff_annotation']
            assert math.isfinite(record['target_score']) and math.isfinite(record['cutoff_score'])
            expected = max(0., R.CONFIG['search_retention']['margin_log_score_units']+record['cutoff_score']-record['target_score']) if annotation['label'] == 0 else 0.
            assert abs(expected-record['margin_loss']) <= 1e-3
            assert disposition == ('known_nonmatch_margin' if annotation['label'] == 0 else 'primary_or_unresolved_zero_margin')
            if annotation['label'] != 0:
                assert record['margin_loss'] == 0.
            qualified_pairs += int(annotation['label'] == 0); terms.append(record['margin_loss'])
        expected_margin = sum(terms)/len(assigned) if assigned else 0.
        assert abs(expected_margin-row['search_retention_margin']) <= 1e-5*(1+abs(expected_margin))
        assert set(row['losses']) == {'local', 'existence_bce', 'boundary_ce', 'termination_bce'}
        assert all(math.isfinite(v) for v in (*row['losses'].values(), row['weighted_total'], row['weighted_original_total'],
                   row['gradient_norm_before_clip'], row['search_retention_margin']))
        weighted_original = row['importance_weight']*sum(row['losses'].values())
        expected_total = weighted_original+row['importance_weight']*R.CONFIG['search_retention']['coefficient']*expected_margin
        assert abs(weighted_original-row['weighted_original_total']) <= 1e-5*(1+abs(weighted_original))
        assert abs(expected_total-row['weighted_total']) <= 1e-5*(1+abs(expected_total))
        margin_steps[species] += int(row['search_retention_margin'] > 0)
        counts[epoch][species].update({'query_'+k: v for k, v in row['query_label_counts'].items()})
        counts[epoch][species].update(row['protection_reasons'])
        losses[epoch].update(dict(row['losses'], search_retention_margin=row['search_retention_margin'], weighted_total=row['weighted_total']))
    assert dict(margin_steps) == fitted['positive_search_margin_steps_by_species'] and all(margin_steps[s] > 0 for s in R.A.ALLOW)
    original = json.loads((R.A.STAGE_DIRS['fit']/'summary.json').read_text()); R.A.check_last_identity(original)
    assert len(fitted['epochs']) == len(original['epochs']) == 3
    for i in range(3):
        actual_supervision = {s: dict(c) for s, c in counts[i].items()}
        assert actual_supervision == fitted['epochs'][i]['supervision_by_species'] == original['epochs'][i]['supervision_by_species']
        assert fitted['epochs'][i]['steps'] == 1536
        assert {k: v/1536 for k, v in losses[i].items()} == fitted['epochs'][i]['mean_losses']
    result = dict(completed=True, arm='A30', source_jobs=dict(fit=fit_job, readout=readout_job), windows=1491,
        actual_4608_draw_weight_and_A3_supervision_identity_verified=True,
        actual_current_search_saved_negative_qualification_and_margin_recomputed=True,
        actual_new_and_original_natural_chain_labels_recomputed=True,
        actual_saved_pool_identity_and_five_fixed_gate_checks_verified=True,
        last_model_structure_and_all_tensors_finite=True, trainable_parameters=649457,
        reference_injected_candidates=0, candidate_window_pairs=candidate_pairs,
        safe_known_nonmatch_pairs=qualified_pairs, positive_query_dispositions=dict(dispositions),
        natural_counts_per_species=actual, original_frozen_TRAIN_generator_gate=checks,
        TRAIN_generator_supported=all(checks.values()), optimizer_steps_this_run=0,
        new_head_or_backbone_forward=False, new_candidate_search=False, DEV_forward_used=False, test_setaria_used=False,
        limitation='Saved loss/label/coverage reconciliation, not a CPU search replay, selected gene F1 or independent generalization.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-job', required=True); parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.fit_job, args.readout_job)
