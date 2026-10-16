#!/usr/bin/env python3
"""Reconcile saved A16 identities and probability decompositions, without forward."""
import json
import math
import torch
import analyze_legal_support_train as D

OUT = D.R.ROOT/'outputs/M28-A16-LEGAL-SUPPORT-RESULT-CHECK'


def main():
    assert (D.OUT/'STATUS').read_text().strip() == 'COMPLETED'
    assert (D.OUT/'JOBID').read_text().strip() == '13367351'
    report = json.loads((D.OUT/'summary.json').read_text())
    assert report['completed'] and report['protocol'] == D.PROTOCOL and report['generator_checkpoint'] == D.C.GENERATOR
    assert report['optimizer_steps_this_run'] == 0 and not report['checkpoint_written']
    assert not report['new_candidate_search'] and not report['new_backbone_forward']
    assert not report['DEV_forward_used'] and not report['test_setaria_used']
    prior = D.R.rows(D.SOURCE/'query_profiles.jsonl')
    source = {(r['window_id'], r['reference_id'], r['query_index']): r for r in prior}
    profiles = D.R.rows(D.OUT/'query_profiles.jsonl')
    keys = [(r['window_id'], r['reference_id'], r['query_index']) for r in profiles]
    assert len(profiles) == len(set(keys)) == len(source) == 432 and set(keys) == set(source)
    maximum_delta = 0.; decisions = 0
    for row, key in zip(profiles, keys):
        old = source[key]; assert row['reference_conditioned'] and len(row['events']) == len(old['events'])
        raw_total = 0.; conditional_total = 0.
        for e, p in zip(row['events'], old['events']):
            assert all(math.isfinite(v) for v in e.values() if isinstance(v, float))
            assert all(e[k] <= 1e-5 for k in ('entry_log_mass', 'end_boundary_log_mass', 'end_joint_log_mass'))
            assert e['end_joint_log_mass'] <= e['end_boundary_log_mass']+1e-5
            assert all(e[k] == p[k] for k in ('exon', 'start', 'end', 'end_kind',
                       'entry_legal_rank', 'end_same_kind_rank', 'carry_before', 'carry_after'))
            maximum_delta = max(maximum_delta, *(abs(e[k]-p[k]) for k in ('entry_logp', 'end_logp', 'eos_logp')))
            assert abs(e['entry_conditional_target_logp']-(e['entry_logp']-e['entry_log_mass'])) < 1e-5
            assert abs(e['end_joint_conditional_target_logp']-(e['end_logp']+e['eos_logp']-e['end_joint_log_mass'])) < 1e-5
            assert e['decoder_vs_teacher_lower_target_logp_change'] == e['decoder_entry_target_logp']-e['teacher_entry_target_logp']
            raw_total += e['entry_logp']+e['end_logp']+e['eos_logp']
            conditional_total += e['entry_conditional_target_logp']+e['end_joint_conditional_target_logp']
            decisions += 1
        assert raw_total == row['raw_path_log_score'] and conditional_total == row['hypothetical_conditioned_path_log_score']
        maximum_delta = max(maximum_delta, abs(raw_total-old['path_log_score']))
        assert row['local_width4_necessary_conditions_passed'] == old['local_width4_necessary_conditions_passed']
    assert maximum_delta <= 1e-3 and maximum_delta == report['raw_score_max_absolute_difference']
    refs = D.R.rows(D.OUT/'references.jsonl'); original = D.R.rows(D.SOURCE/'references.jsonl')
    assert len(refs) == len(original) == 27
    boundary_means = []
    for row, old in zip(refs, original):
        assert {k: v for k, v in row.items() if k != 'support'} == old
        local = [p for p in profiles if (p['window_id'], p['reference_id']) == (row['window_id'], row['reference_id'])]
        assert len(local) == 16 and row['support'] == D.reference_means(local)
        events = [e for p in local for e in p['events']]
        boundary_means.append(sum(math.exp(e['end_boundary_log_mass']) for e in events)/len(events))
    assert report['reference_balanced_pooled'] == D.aggregate(refs)
    assert report['per_species'] == {s: D.aggregate([r for r in refs if r['species'] == s]) for s in D.R.ALLOW}
    assert report['by_exon_class'] == {k: D.aggregate([r for r in refs if r['exon_class'] == k]) for k in ('one', 'two_to_four', 'five_or_more')}
    assert report['by_original_natural_beam_hit'] == {str(h): D.aggregate([r for r in refs if r['beam_free_hit'] == h]) for h in (False, True)}
    assert report['query_profiles'] == 432 and report['reference_window_pairs'] == 27 and report['exon_decisions'] == decisions
    result = dict(completed=True, source_job='13367351', windows=8, query_profiles=432,
        reference_window_pairs=27, exon_decisions=decisions,
        actual_source_identity_local_ranks_carry_and_raw_scores_verified=True,
        actual_conditional_probability_decomposition_and_aggregates_verified=True,
        raw_score_max_absolute_difference=maximum_delta,
        reference_balanced_pooled=report['reference_balanced_pooled'],
        per_species=report['per_species'], by_exon_class=report['by_exon_class'],
        by_original_natural_beam_hit=report['by_original_natural_beam_hit'],
        end_boundary_probability_mass_reference_balanced=sum(boundary_means)/len(boundary_means),
        optimizer_steps_this_run=0, new_head_or_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False,
        limitation='Saved TRAIN reference-conditioned probe only; no new natural predictions, model benefit or causal attribution.')
    with (OUT/'summary.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2); main()
