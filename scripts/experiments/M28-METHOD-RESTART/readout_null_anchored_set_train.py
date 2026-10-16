#!/usr/bin/env python3
"""Last A39 model on the complete original label-free natural TRAIN pool."""
import argparse
from collections import Counter
import json
import time
import torch
import fit_null_anchored_set as F
import readout_logistic_context_train as L

ROOT, CONFIG = F.ROOT, F.CONFIG
B, R = F.A.B, F.A.B.R
FIT_JOB = None
PROTOCOL = F.PROTOCOL
OUT = ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-RESULT-CHECK'
CHECKS = ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-READOUT-CHECKS-R1'
COMPONENTS = ('additive', 'nonadditive', 'exon_additive', 'junction_additive')


def validate_fit(value, state=False):
    if (value['arm'] != 'A39' or value['configuration'] != CONFIG or value['protocol'] != F.PROTOCOL or
            value['step' if state else 'optimizer_steps'] != 4608 or
            value['revision' if state else 'feature_revision'] != F.A.R.C.I.F.REVISION):
        raise ValueError('Only the fixed complete A39 last004608 model is allowed')
    if not state and (value['completed'] is not True or value['preflight_job'] != F.PREFLIGHT_JOB or
            value['trainable_parameters'] != 612324 or value['initialization'] != CONFIG['initialization'] or
            value['final_checkpoint'] != str((F.OUT/'step_004608.pt').relative_to(ROOT)) or
            value['reference_injected_training_candidates'] != 0 or
            any(value[k] is not False for k in ('supervised_checkpoint_read','generator_forward',
                'generator_updated','new_candidate_search','new_backbone_forward','new_R4_labelling',
                'DEV_used_for_selection','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual complete unchanged natural-only A39 fit required')


def gate(counts):
    checks = {}
    for species, values in counts.items():
        old, new = values['A31_A3'], values['A39']
        checks[species+'_natural_selected_correct_strictly_increased'] = new['selected_correct'] > old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = (
            new['selected_known_nonmatch'] <= old['selected_known_nonmatch'])
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A39']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def paired_readout(saved, previous, labelled, context, logistic):
    if not len(saved) == len(previous) == len(labelled) == len(context) == len(logistic):
        raise ValueError('Complete paired original windows required')
    totals = {s: {a: Counter() for a in ('A31_A3','A37','A38','A39')} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    rows = []
    for new, old, lab, prior_context, prior_logistic in zip(saved, previous, labelled, context, logistic):
        if (new['window'] != old['window'] or new['window'] != lab['window'] or
                new['reference_used_in_scoring'] or prior_context['reference_used_in_scoring'] or
                prior_logistic['reference_used_in_scoring'] or
                new['window'] != prior_logistic['window'] or new['chains'] != prior_logistic['chains'] or
                new['window'] != prior_context['window'] or new['chains'] != prior_context['chains'] or
                new['chains'] != old['variants']['A31_A3']['chains']):
            raise ValueError('Unchanged ordered reference-free original natural geometry required')
        labels = lab['labels']['A31_A3']
        if any(y not in (-1,0,1) for y in labels):
            raise ValueError('Original protected labels required')
        counts = dict(A31_A3=R.selected_counts(old['variants']['A31_A3'], labels),
                      A37=R.selected_counts(prior_context, labels), A38=R.selected_counts(prior_logistic, labels),
                      A39=R.selected_counts(new, labels))
        if counts['A31_A3'] != lab['counts']['A31_A3']:
            raise ValueError('Original stronger-baseline selection differs')
        species = new['window']['species']
        for arm in totals[species]:
            totals[species][arm].update(counts[arm])
        references[species].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
        rows.append(dict(window=new['window'], labels=labels, counts=counts,
            eligible_references=lab['eligible_references'], eligible_five_or_more=lab['eligible_five_or_more']))
    return rows, {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}, {
        s: dict(c) for s, c in references.items()}


def readiness():
    F.A.closed_a38()
    draws, free, previous, index, old = B.population()
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    context = R.C.R.rows(L.F.A.D.OUT/'natural_scores.jsonl')
    context_checked = R.completed(L.F.A.D.CHECK,'13410742')
    logistic = R.C.R.rows(L.OUT/'natural_scores.jsonl')
    logistic_checked = R.completed(L.CHECK, '13412702')
    baseline = [dict(window=r['window'],reference_used_in_scoring=False,
                     **r['variants']['A31_A3']) for r in previous]
    rows, counts, references = paired_readout(baseline,previous,labelled,context,logistic)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A39'] != counts[s]['A31_A3'] or
                counts[s]['A37'] != context_checked['selected_counts_per_species'][s]['A37'] or
                counts[s]['A38'] != logistic_checked['selected_counts_per_species'][s]['A38'] for s in counts)):
        raise ValueError('Original complete four-arm metadata pairing differs')
    R.C.R.write_summary(CHECKS,dict(completed=True,experiment='M28-A39-NULL-ANCHORED-SET-READOUT-CHECKS',
        protocol=PROTOCOL,configuration=CONFIG,windows=1491,original_draws=len(draws),
        natural_candidates=sum(len(v['chains']) for v in free),feature_index_records=len(index),
        original_four_arm_identity_labels_DP_baseline_and_denominators_reconciled=True,
        new_fitted_checkpoint_opened=False,feature_tensors_opened=False,model_or_backbone_forward=False,
        optimizer_steps_this_run=0,new_R4_labelling=False,new_candidate_search=False,DEV_used=False,
        test_setaria_used=False,A33_outputs_opened=False,
        limitation='Original metadata and synthetic pairing checks only; no A39 performance.'))


@torch.inference_mode()
def main(checks_job):
    F.A.require_jobs([FIT_JOB, F.PREFLIGHT_JOB,str(checks_job)])
    prepared = R.completed(CHECKS,checks_job)
    if (prepared['protocol'] != PROTOCOL or prepared['configuration'] != CONFIG or
            prepared['original_four_arm_identity_labels_DP_baseline_and_denominators_reconciled'] is not True):
        raise ValueError('Actual unchanged four-arm pairing readiness required')
    F.ready()
    fitted = R.completed(F.OUT, FIT_JOB); validate_fit(fitted)
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    model = F.A.new_model('cuda').eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True); del state
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    draws, free, previous, index, old = B.population()
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal, pool in enumerate(free):
            w = pool['window']; cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            score = model(features, dna, [w['valid_bases']], pool['chains'])
            if any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
                raise ValueError('Complete finite FP32 joint-head scores required')
            values = score['logits'].cpu().tolist()
            row = dict(window=w, chains=pool['chains'], decoder_sources=pool['decoder_sources'], scores=values,
                selected=R.select_nonoverlapping(pool['chains'], values),
                **{k: score[k].cpu().tolist() for k in COMPONENTS}, null=float(score['null']),
                reference_used_in_scoring=False)
            handle.write(json.dumps(row, allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0 or ordinal+1 == len(free):
                handle.flush(); print(json.dumps(dict(natural_windows=ordinal+1, total=1491)), flush=True)
            del cached, features, dna, score
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen readout changed weights or received gradients')
    # Original saved labels join only after every new natural score is closed.
    saved = R.C.R.rows(OUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    context = R.C.R.rows(L.F.A.D.OUT/'natural_scores.jsonl')
    context_checked = R.completed(L.F.A.D.CHECK,'13410742')
    logistic = R.C.R.rows(L.OUT/'natural_scores.jsonl')
    logistic_checked = R.completed(L.CHECK, '13412702')
    rows, counts, references = paired_readout(saved, previous, labelled, context, logistic)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A37'] != context_checked['selected_counts_per_species'][s]['A37'] or
                counts[s]['A38'] != logistic_checked['selected_counts_per_species'][s]['A38'] for s in counts)):
        raise ValueError('Original full denominators and stronger baseline differ')
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    checks = gate(counts)
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A39-NULL-ANCHORED-SET-TRAIN-READOUT',
        protocol=PROTOCOL, fit_protocol=F.PROTOCOL, configuration=CONFIG, source_fit_job=FIT_JOB,
        source_preflight_job=F.PREFLIGHT_JOB,source_checks_job=str(checks_job),
        source_context_comparator_jobs=dict(fit='13410061',readout='13410739',CPU='13410742'),
        source_logistic_comparator_jobs=dict(fit='13412630',readout='13412701',CPU='13412702'), source_baseline_jobs=dict(readout='13398521', CPU='13398524'),
        checkpoint=fitted['final_checkpoint'], windows=1491, original_draws=1536, natural_candidates=16321,
        reference_injected_predictions=0, natural_score_file_closed_before_saved_label_arrays_open=True,
        original_natural_geometry_order_sources_unchanged=True, original_saved_labels_reused=True,
        label_catalog_opened=False, selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0, new_candidate_search=False, new_backbone_forward=False,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='In-sample original natural window/chain decisions only; not free generation, '
            'independent gene F1, causal mechanism or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('checks','readout'),required=True)
    parser.add_argument('--fit-job')
    parser.add_argument('--checks-job')
    args = parser.parse_args(); FIT_JOB = args.fit_job
    torch.set_num_threads(2)
    if args.mode == 'checks':
        readiness()
    else:
        if not args.fit_job or not args.checks_job:
            parser.error('--fit-job and --checks-job are required for actual readout')
        main(args.checks_job)
