#!/usr/bin/env python3
"""A31 CPU-only saved-pool, R4-label, selection and five-gate reconciliation."""
import argparse
import json
import torch
import run_local_path_caller as R


def main(source_job):
    report = R.completed(R.READOUT, source_job)
    if (report['configuration'] != R.CONFIG or report['protocol'] != R.PROTOCOL or
            report['windows'] != 1491 or report['optimizer_steps_this_run'] != 0 or
            not report['natural_score_file_closed_before_catalog_open'] or
            not report['original_greedy_source_replay_verified'] or
            not report['model_weights_unchanged'] or not report['all_parameters_frozen_and_gradients_absent'] or
            report['new_candidate_search'] or report['new_backbone_forward'] or
            report['DEV_forward_used'] or report['test_setaria_used'] or report['reference_injected_candidates'] != 0):
        raise ValueError('Matching completed frozen TRAIN-only A31 readout required')
    draws, ids, index, pools = R.sources()
    loaded, _, windows, catalogs = R.Q.load_train_inputs()
    saved = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    if loaded != draws or [r['window']['id'] for r in saved] != ids or [r['window']['id'] for r in labelled] != ids:
        raise ValueError('Complete ordered original TRAIN records required')
    sequences = {}
    for ident in ids:
        cached = torch.load(R.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        R.Q.validate_cached(cached, windows[ident]); sequences[ident] = cached['sequence']
        del cached
    auxiliary = {r['window']['id']: r for r in R.C.R.rows(R.D.OUT/'windows.jsonl')}
    actual, counts, references = R.label_readout(saved, pools, windows, catalogs, sequences, auxiliary)
    if (actual != labelled or counts != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or
            R.proxies(counts, references) != report['known_label_window_pair_proxies'] or
            sum(v['eligible'] for v in references.values()) != 2916):
        raise ValueError('All saved labels, selections, reference denominators and proxy metrics must reconcile')
    checks = R.gate(counts, references)
    if checks != report['predeclared_TRAIN_gate'] or {
            a: all(v.values()) for a, v in checks.items()} != report['TRAIN_preparation_supported']:
        raise ValueError('Original five-gate integer arithmetic differs')
    for a, path in R.CONFIG['generators'].items():
        state = torch.load(R.ROOT/path, map_location='cpu', weights_only=True)
        (R.C.A.check_last_identity if a == 'A3' else R.R23.validate_fit)(state, state=True)
        if (not all(bool(torch.isfinite(v).all()) for v in state['model'].values()) or
                sum(v.numel() for v in state['model'].values()) != {'A3': 649457, 'A23': 649617}[a]):
            raise ValueError('Exact finite original last-checkpoint structure required')
        del state
    R.write(R.CHECK, dict(completed=True, source_jobs=dict(readout=str(source_job)), windows=1491,
        actual_saved_natural_pool_and_greedy_subset_identity_verified=True,
        actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed=True,
        actual_five_fixed_gate_integer_arithmetic_verified=True,
        original_last_checkpoint_tensors_finite_and_parameter_counts_verified=True,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported={a: all(v.values()) for a, v in checks.items()},
        new_head_or_backbone_forward=False,
        limitation='CPU reconciles saved emissions/selection and R4 labels; no independent emission forward, gene F1 or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.source_job)
