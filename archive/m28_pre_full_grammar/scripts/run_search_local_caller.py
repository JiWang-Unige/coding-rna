#!/usr/bin/env python3
"""A32: gated frozen A30 pools/model with the unchanged A31 local caller."""
import argparse
from collections import Counter
import json
import os
import subprocess
import time
import torch
import run_local_path_caller as R
import reconcile_query_search_retention_fit as P

READOUT = R.ROOT/'outputs/M28-A32-SEARCH-LOCAL-CALLER-TRAIN-READOUT'
CHECK = R.ROOT/'outputs/M28-A32-SEARCH-LOCAL-CALLER-RESULT-CHECK'
PROTOCOL = 'reports/M28-METHOD-RESTART/a32_search_local_caller_protocol.md'
SOURCE_JOBS = dict(A30_reconciliation='13399750', A30_readout='13399752', A30_CPU='13399755',
                   A31_readout='13398521', A31_CPU='13398524')
CONFIG = dict(generator_checkpoint=str((P.K.SOURCE/'step_004608.pt').relative_to(R.ROOT)),
    generator_optimizer_steps=4608, generator_parameters=649457,
    score=R.CONFIG['score'], all_role_weights=1., null=0., model_precision='FP32',
    score_accumulation='FP64', primary_baseline='A31_A3', threshold=0., threshold_scan=False,
    pool='reuse_actual_closed_A30_greedy_plus_beam4_union', beam_width=4,
    completed_limit=4, max_exons=128, unique_TRAIN_windows=1491,
    new_parameters=0, optimizer_steps=0)
ARMS = ('A31_A3', 'A32_A30')


def require_roots(jobs):
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    receipts = {}
    for line in result.stdout.splitlines():
        fields = [v.strip() for v in line.split('|')]
        if len(fields) >= 3 and fields[0] in jobs:
            receipts[fields[0]] = fields[1:3]
    if any(receipts.get(job) != ['COMPLETED', '0:0'] for job in jobs):
        raise ValueError('Actual completed source roots required before A32 execution')


def require_generator_gate(report, checked):
    if (not checked['actual_saved_pool_identity_and_five_fixed_gate_checks_verified'] or
            not checked['actual_new_and_original_natural_chain_labels_recomputed'] or
            checked['natural_counts_per_species'] != report['natural_counts_per_species'] or
            set(checked['natural_counts_per_species']) != set(R.C.R.ALLOW)):
        raise ValueError('Actual original A30 pool/label/gate reconciliation required')
    checks = P.R.gate(checked['natural_counts_per_species'])
    if (len(checks) != 5 or checks != report['predeclared_TRAIN_generator_gate'] or
            checks != checked['original_frozen_TRAIN_generator_gate'] or
            all(checks.values()) != report['TRAIN_generator_supported'] or
            all(checks.values()) != checked['TRAIN_generator_supported']):
        raise ValueError('Original A30 five-gate accounting differs')
    if not all(checks.values()):
        raise ValueError('Original A30 generator gate failed; A32 not executed')


def sources():
    require_roots(list(SOURCE_JOBS.values()))
    P.require_reconciled(SOURCE_JOBS['A30_reconciliation'])
    fitted = R.completed(P.K.SOURCE, P.K.SOURCE_JOB); P.R.validate_fit(fitted)
    parent = R.completed(P.READOUT, SOURCE_JOBS['A30_readout'])
    checked = R.completed(P.CHECK, SOURCE_JOBS['A30_CPU'])
    if (parent['configuration'] != P.R.CONFIG or parent['protocol'] != P.R.PROTOCOL or
            parent['checkpoint'] != CONFIG['generator_checkpoint'] or parent['windows'] != 1491 or
            checked['source_jobs'] != dict(fit=P.K.SOURCE_JOB, readout=SOURCE_JOBS['A30_readout']) or
            checked['source_fit_root_state'] != 'TIMEOUT' or
            checked['source_fit_artifact_reconciliation_job'] != SOURCE_JOBS['A30_reconciliation'] or
            checked['actual_TRAIN_optimizer_updates_this_recovery'] != 0):
        raise ValueError('Matching original A30 final model and readout-only recovery required')
    require_generator_gate(parent, checked)
    baseline = R.completed(R.READOUT, SOURCE_JOBS['A31_readout'])
    baseline_check = R.completed(R.CHECK, SOURCE_JOBS['A31_CPU'])
    if (baseline['configuration'] != R.CONFIG or baseline_check['configuration'] != R.CONFIG or
            baseline['protocol'] != R.PROTOCOL or baseline_check['protocol'] != R.PROTOCOL or
            baseline_check['source_jobs'] != dict(readout=SOURCE_JOBS['A31_readout']) or
            not baseline_check['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'] or
            not baseline_check['TRAIN_preparation_supported']['A31_A3'] or
            baseline['selected_counts_per_species'] != baseline_check['selected_counts_per_species'] or
            baseline['eligible_reference_pairs'] != baseline_check['eligible_reference_pairs']):
        raise ValueError('Original complete A31 A3 caller baseline required')
    return fitted, parent, baseline


def gate(counts, references):
    checks = {}
    for species, values in counts.items():
        new, old = values['A32_A30'], values['A31_A3']
        n, o = new['selected_correct'], old['selected_correct']
        nd = references[species]['eligible']+n+new['selected_known_nonmatch']
        od = references[species]['eligible']+o+old['selected_known_nonmatch']
        checks[species+'_proxy_F1_strictly_increased'] = n*od > o*nd
        checks[species+'_proxy_precision_not_decreased'] = R.precision_compare(new, old)
    checks['pooled_five_or_more_selected_correct_strictly_increased'] = (
        sum(v['A32_A30']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def paired_counts(scored, raw, labelled, baseline, baseline_labelled):
    for other in (raw, labelled, baseline, baseline_labelled):
        if any(scored['window'][k] != other['window'][k] for k in R.C.R.GEOMETRY):
            raise ValueError('Original paired TRAIN geometry differs')
    if (scored['reference_used_in_scoring'] or raw['reference_used_in_pool_generation'] or
            scored['chains'] != raw['chains'] or
            scored['chains'] != [c['chain'] for c in labelled['candidates']]):
        raise ValueError('Do not change closed A30 pool geometry or label order')
    labels = [c['label'] for c in labelled['candidates']]
    new = R.selected_counts(scored, labels)
    old = R.selected_counts(baseline['variants']['A31_A3'], baseline_labelled['labels']['A31_A3'])
    if old != baseline_labelled['counts']['A31_A3']:
        raise ValueError('Original A31 A3 per-window counts differ')
    return dict(A31_A3=old, A32_A30=new), labels


def tabulate(saved, parent, baseline):
    raw = R.C.R.rows(P.READOUT/'free_pool.jsonl')
    labelled = R.C.R.rows(P.READOUT/'windows.jsonl')
    old = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    old_labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    ids = list(dict.fromkeys(d['window_id'] for d in R.C.R.rows(R.Q.DRAWS)))
    if len(ids) != 1491 or any([r['window']['id'] for r in rows] != ids
                              for rows in (saved, raw, labelled, old, old_labelled)):
        raise ValueError('Complete original ordered 1491-window population required')
    counts = {s: {a: Counter() for a in ARMS} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}; records = []
    for value, pool, labels, prior, prior_labels in zip(saved, raw, labelled, old, old_labelled):
        paired, new_labels = paired_counts(value, pool, labels, prior, prior_labels)
        species = value['window']['species']
        for arm in ARMS:
            counts[species][arm].update(paired[arm])
        references[species].update(eligible=prior_labels['eligible_references'],
                                   five_or_more=prior_labels['eligible_five_or_more'])
        records.append(dict(window=value['window'], labels=new_labels, counts=paired,
            eligible_references=prior_labels['eligible_references'],
            eligible_five_or_more=prior_labels['eligible_five_or_more']))
    counts = {s: {a: dict(c) for a, c in v.items()} for s, v in counts.items()}
    references = {s: dict(c) for s, c in references.items()}
    if references != baseline['eligible_reference_pairs'] or sum(v['eligible'] for v in references.values()) != 2916:
        raise ValueError('All original eligible reference pairs must remain in the denominator')
    for species, values in counts.items():
        if values['A31_A3'] != baseline['selected_counts_per_species'][species]['A31_A3']:
            raise ValueError('Complete original A31 A3 baseline must reproduce exactly')
        expected = parent['natural_counts_per_species'][species]['A30']
        if tuple(values['A32_A30'][k] for k in ('candidates', 'pool_correct', 'pool_known_nonmatch', 'pool_unknown')) != \
                tuple(expected[k] for k in ('free_candidates', 'correct', 'known_nonmatch', 'unknown')):
            raise ValueError('Actual closed A30 natural-pool totals differ')
    return records, counts, references


def write(folder, **result):
    R.C.R.write_summary(folder, dict(completed=True, arm='A32_A30', configuration=CONFIG, protocol=PROTOCOL,
        source_jobs=SOURCE_JOBS, optimizer_steps_this_run=0, new_parameters=0, new_candidate_search=False,
        new_backbone_forward=False, new_R4_labelling=False, reference_injected_candidates=0,
        DEV_used_for_selection=False, test_setaria_used=False, **result))


def load_model(fitted):
    state = torch.load(R.ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    P.R.validate_fit(state, state=True)
    model = P.R.new_model().eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True)
    if sum(p.numel() for p in model.parameters()) != CONFIG['generator_parameters']:
        raise ValueError('Unchanged A30/A3 capacity required')
    return model.to('cuda')


@torch.inference_mode()
def readout():
    # Fail the original source gate before model weights, feature cache or pool are opened.
    fitted, parent, baseline = sources()
    pools = R.C.R.rows(P.READOUT/'free_pool.jsonl')
    ids = list(dict.fromkeys(d['window_id'] for d in R.C.R.rows(R.Q.DRAWS)))
    if len(ids) != 1491 or [p['window']['id'] for p in pools] != ids:
        raise ValueError('Original ordered 1491 TRAIN windows required')
    index = {r['window_id']: r for r in R.Q.train_records(R.Q.FEATURES/'index.jsonl')
             if r['window_id'] in set(ids)}
    if set(index) != set(ids):
        raise ValueError('Original TRAIN feature population required')
    model = load_model(fitted)
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    READOUT.mkdir(exist_ok=False)
    (READOUT/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (READOUT/'STATUS').write_text('RUNNING\n')
    began = time.perf_counter(); saved = []; torch.cuda.reset_peak_memory_stats()
    # Scores use only the closed natural pool and the frozen model, never saved R4 labels.
    with (READOUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal, pool in enumerate(pools):
            w = pool['window']
            if (pool['reference_used_in_pool_generation'] or w['split'] != 'train' or
                    w['seqid'] not in R.C.R.ALLOW[w['species']]):
                raise ValueError('Reference-free allowed TRAIN geometry required')
            cached = torch.load(R.ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
            R.Q.validate_cached(cached, w)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]
            logits = model.segmentation(hidden)
            if logits.dtype != torch.float32:
                raise ValueError('Unchanged frozen FP32 forward required')
            scores = R.local_path_scores(logits, pool['chains'], w['valid_bases']).cpu().tolist()
            value = dict(window=w, chains=pool['chains'], scores=scores,
                selected=R.select_nonoverlapping(pool['chains'], scores), reference_used_in_scoring=False)
            handle.write(json.dumps(value, allow_nan=False)+'\n'); saved.append(value)
            if (ordinal+1) % 64 == 0:
                handle.flush(); print(json.dumps(dict(phase='natural_scores', windows=ordinal+1, total=1491)), flush=True)
            del cached, features, dna, hidden, logits
    torch.cuda.synchronize(); scoring_seconds = time.perf_counter()-began
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen model changed or received a gradient')
    labelled, counts, references = tabulate(saved, parent, baseline)
    with (READOUT/'windows.jsonl').open('x') as handle:
        for value in labelled:
            handle.write(json.dumps(value, allow_nan=False)+'\n')
    checks = gate(counts, references)
    write(READOUT, windows=1491, generator_gate_passed_before_weights_or_cache_open=True,
        natural_score_file_closed_before_saved_per_window_labels_open=True,
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        known_label_window_pair_proxies=R.proxies(counts, references), predeclared_TRAIN_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), scoring_seconds=scoring_seconds,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One-seed matched-capacity/exposure TRAIN composite, not independent causality, '
            'full gene F1, native competitiveness, automatic DEV release or generalization.')
    (READOUT/'STATUS').write_text('COMPLETED\n')


def check(source_job):
    _, parent, baseline = sources(); require_roots([source_job])
    report = R.completed(READOUT, source_job)
    if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
            not report['model_weights_unchanged'] or
            not report['all_parameters_frozen_and_gradients_absent'] or
            not report['generator_gate_passed_before_weights_or_cache_open'] or
            not report['natural_score_file_closed_before_saved_per_window_labels_open']):
        raise ValueError('Matching frozen complete A32 natural score file required')
    saved = R.C.R.rows(READOUT/'natural_scores.jsonl')
    labelled, counts, references = tabulate(saved, parent, baseline)
    checks = gate(counts, references)
    if (labelled != R.C.R.rows(READOUT/'windows.jsonl') or
            counts != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or
            R.proxies(counts, references) != report['known_label_window_pair_proxies'] or
            checks != report['predeclared_TRAIN_gate'] or
            all(checks.values()) != report['TRAIN_preparation_supported']):
        raise ValueError('Actual A32 saved selection, denominators or five-gate accounting differs')
    CHECK.mkdir(exist_ok=False)
    (CHECK/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (CHECK/'STATUS').write_text('RUNNING\n')
    write(CHECK, windows=1491, readout_job=source_job,
        actual_saved_pool_and_upstream_R4_labels_selection_denominators_reconciled=True,
        actual_five_fixed_gate_integer_arithmetic_verified=True,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False,
        limitation='Uses independently validated upstream R4 labels; no independent emission replay, '
            'new R4 labelling, raw-genome annotation or biological generalization.')
    (CHECK/'STATUS').write_text('COMPLETED\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('readout', 'check'), required=True)
    parser.add_argument('--source-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'readout':
        readout()
    elif not args.source_job:
        raise ValueError('Actual A32 readout job ID required')
    else:
        check(args.source_job)
