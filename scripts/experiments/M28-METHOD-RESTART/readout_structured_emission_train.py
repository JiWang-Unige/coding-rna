#!/usr/bin/env python3
"""A36 last-model readout on the original natural pool, never injected chains."""
from collections import Counter
import argparse
import json
import time
import torch
import fit_structured_emission as F
import run_local_path_caller as R
from src.m28.local_path_score import local_path_scores

ROOT, CONFIG = F.ROOT, F.CONFIG
OUT = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-RESULT-CHECK'
CHECKS = ROOT/'outputs/M28-A36-STRUCTURED-EMISSION-READOUT-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a36_structured_emission_readout_protocol.md'
FIT_JOB = '13406144'


def validate_fit(value, state=False):
    if (value['arm'] != 'A36' or value['configuration'] != CONFIG or value['protocol'] != F.PROTOCOL or
            value['step' if state else 'optimizer_steps'] != 4608 or
            value['revision' if state else 'feature_revision'] != F.A.C.I.F.REVISION):
        raise ValueError('Only fixed A36 last004608 is permitted')
    if not state and (not value['completed'] or value['preflight_job'] != F.PREFLIGHT_JOB or
            value['final_checkpoint'] != str((F.OUT/'step_004608.pt').relative_to(ROOT)) or
            value['trainable_parameters'] != 366191 or value['DEV_used_for_selection'] or
            value['test_setaria_used'] or value['A33_outputs_opened']):
        raise ValueError('Actual complete unchanged TRAIN-only A36 fit required')


def population():
    # No TRAIN catalog, saved label arrays or fitted new model is opened here.
    F.A.require_jobs(['13398521', '13398524'])
    old = R.completed(R.READOUT, 13398521); checked = R.completed(R.CHECK, 13398524)
    if (old['configuration'] != R.CONFIG or checked['configuration'] != R.CONFIG or
            old['protocol'] != R.PROTOCOL or checked['protocol'] != R.PROTOCOL or
            old['windows'] != 1491 or checked['windows'] != 1491 or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] or
            not checked['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'] or
            old['selected_counts_per_species'] != checked['selected_counts_per_species']):
        raise ValueError('Actual independently reconciled original A31_A3 baseline required')
    draws = R.C.R.rows(R.Q.DRAWS); ids = list(dict.fromkeys(d['window_id'] for d in draws))
    free = R.C.R.rows(R.C.FIT/'free_pool.jsonl')
    prior = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    if (len(draws) != 1536 or len(ids) != 1491 or
            [r['window']['id'] for r in free] != ids or [r['window']['id'] for r in prior] != ids or
            sum(len(r['chains']) for r in free) != 16321):
        raise ValueError('Complete ordered original natural pool required')
    for pool, native in zip(free, prior):
        w = pool['window']
        if (w['split'] != 'train' or w['seqid'] not in R.C.R.ALLOW[w['species']] or
                pool['reference_used_in_pool_generation'] or native['reference_used_in_scoring'] or
                w != native['window'] or pool['chains'] != native['variants']['A31_A3']['chains'] or
                len(pool['decoder_sources']) != len(pool['chains'])):
            raise ValueError('Unchanged reference-free A3 natural geometry and provenance required')
    idset = set(ids)
    index = {r['window_id']: r for r in R.Q.train_records(R.Q.FEATURES/'index.jsonl')
             if r['window_id'] in idset}
    if set(index) != idset:
        raise ValueError('Original TRAIN feature coverage required')
    return draws, free, prior, index, old


def paired_readout(saved, prior, labelled):
    if len(saved) != len(prior) or len(saved) != len(labelled):
        raise ValueError('Complete paired original windows required')
    totals = {s: {a: Counter() for a in ('A31_A3', 'A36')} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    rows = []
    for new, old, lab in zip(saved, prior, labelled):
        if (new['window'] != old['window'] or new['window'] != lab['window'] or
                new['reference_used_in_scoring'] or new['chains'] != old['variants']['A31_A3']['chains']):
            raise ValueError('No reordered, added or reference-injected prediction geometry')
        labels = lab['labels']['A31_A3']
        if any(y not in (-1, 0, 1) for y in labels):
            raise ValueError('Original protected 1/0/-1 labels required')
        counts = dict(A31_A3=R.selected_counts(old['variants']['A31_A3'], labels),
                      A36=R.selected_counts(new, labels))
        if counts['A31_A3'] != lab['counts']['A31_A3']:
            raise ValueError('Original A31_A3 counts differ')
        species = new['window']['species']
        for arm in totals[species]:
            totals[species][arm].update(counts[arm])
        references[species].update(eligible=lab['eligible_references'],
                                   five_or_more=lab['eligible_five_or_more'])
        rows.append(dict(window=new['window'], labels=labels, counts=counts,
            eligible_references=lab['eligible_references'], eligible_five_or_more=lab['eligible_five_or_more']))
    return rows, {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}, {
        s: dict(c) for s, c in references.items()}


def gate(counts):
    checks = {}
    for species, values in counts.items():
        new, old = values['A36'], values['A31_A3']
        checks[species+'_natural_selected_correct_strictly_increased'] = new['selected_correct'] > old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = (
            new['selected_known_nonmatch'] <= old['selected_known_nonmatch'])
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A36']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def readiness():
    # Reconcile existing baseline labels only, without opening the new fit.
    draws, free, prior, index, old = population()
    labels = R.C.R.rows(R.READOUT/'windows.jsonl')
    baseline = [dict(window=r['window'], reference_used_in_scoring=False,
                     **r['variants']['A31_A3']) for r in prior]
    rows, counts, references = paired_readout(baseline, prior, labels)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A36'] != counts[s]['A31_A3'] for s in counts)):
        raise ValueError('Original complete identity/label join or stronger baseline differs')
    R.C.R.write_summary(CHECKS, dict(completed=True,
        experiment='M28-A36-STRUCTURED-EMISSION-READOUT-CHECKS', protocol=PROTOCOL,
        configuration=CONFIG, windows=len(free), original_draws=len(draws),
        natural_candidates=sum(len(r['chains']) for r in free), feature_index_records=len(index),
        original_identity_labels_stronger_baseline_and_reference_denominators_reconciled=True,
        new_fitted_checkpoint_opened=False, feature_tensors_opened=False,
        model_or_backbone_forward=False, optimizer_steps_this_run=0,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Original saved metadata readiness and synthetic implementation checks only; '
            'no new model scoring or performance result.'))


@torch.inference_mode()
def main(checks_job):
    F.A.require_jobs([str(checks_job)])
    ready = R.completed(CHECKS, checks_job)
    if (ready['protocol'] != PROTOCOL or ready['configuration'] != CONFIG or
            not ready['original_identity_labels_stronger_baseline_and_reference_denominators_reconciled']):
        raise ValueError('Actual matching original metadata readiness required')
    F.A.require_jobs([FIT_JOB]); fitted = R.completed(F.OUT, FIT_JOB); validate_fit(fitted)
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    model = F.A.new_model('cuda').eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True); del state
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    draws, free, prior, index, old = population()
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal, pool in enumerate(free):
            w = pool['window']; ident = w['id']
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            logits = model(features, dna, [w['valid_bases']])['segmentation_logits'][0][:w['valid_bases']]
            if logits.dtype != torch.float32:
                raise ValueError('Original FP32 model precision required')
            scores = local_path_scores(logits, pool['chains']).cpu().tolist()
            selected = R.select_nonoverlapping(pool['chains'], scores)
            path = OUT/'logits'/('%04d.pt' % ordinal)
            with path.open('xb') as file:
                torch.save(dict(window=w, logits=logits.cpu()), file)
            row = dict(window=w, chains=pool['chains'], decoder_sources=pool['decoder_sources'],
                scores=scores, selected=selected, logits_artifact=str(path.relative_to(ROOT)),
                reference_used_in_scoring=False)
            handle.write(json.dumps(row, allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0 or ordinal+1 == len(free):
                handle.flush(); print(json.dumps(dict(natural_windows=ordinal+1, total=1491)), flush=True)
            del cached, features, dna, logits
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen readout altered parameters or received gradients')
    # Join the already independently verified labels only after all scores close.
    saved = R.C.R.rows(OUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    rows, counts, references = paired_readout(saved, prior, labelled)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Original reference denominators or stronger baseline drifted')
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    checks = gate(counts)
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A36-STRUCTURED-EMISSION-TRAIN-READOUT',
        protocol=PROTOCOL, fit_protocol=F.PROTOCOL, configuration=CONFIG, source_fit_job=FIT_JOB,
        source_checks_job=str(checks_job),
        source_baseline_jobs=dict(readout='13398521', CPU='13398524'), checkpoint=fitted['final_checkpoint'],
        windows=1491, original_draws=1536, natural_candidates=16321, reference_injected_predictions=0,
        natural_score_file_closed_before_saved_label_arrays_open=True, original_natural_geometry_order_sources_unchanged=True,
        original_saved_labels_reused=True, label_catalog_opened=False, new_R4_labelling=False,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0, new_candidate_search=False, new_backbone_forward=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='In-sample same-natural-pool window/chain counts only; not free generation, '
            'independent gene F1, causal mechanism, generalization or publication readiness.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('checks', 'readout'), required=True)
    parser.add_argument('--checks-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'checks':
        readiness()
    elif args.checks_job is None:
        parser.error('--checks-job is required for readout')
    else:
        main(args.checks_job)
