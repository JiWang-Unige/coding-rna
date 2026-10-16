#!/usr/bin/env python3
"""Last A37 model on the complete original label-free natural TRAIN pool."""
import argparse
from collections import Counter
import json
import time
import torch
import fit_context_chain as F

ROOT, CONFIG = F.ROOT, F.CONFIG
B, R = F.A.B, F.A.B.R
FIT_JOB = '13410061'
PROTOCOL = 'reports/M28-METHOD-RESTART/a37_context_chain_readout_protocol.md'
OUT = ROOT/'outputs/M28-A37-CONTEXT-CHAIN-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A37-CONTEXT-CHAIN-RESULT-CHECK'
CHECKS = ROOT/'outputs/M28-A37-CONTEXT-CHAIN-READOUT-CHECKS'
COMPONENTS = ('additive', 'nonadditive', 'exon_additive', 'junction_additive')


def validate_fit(value, state=False):
    if (value['arm'] != 'A37' or value['configuration'] != CONFIG or value['protocol'] != F.PROTOCOL or
            value['step' if state else 'optimizer_steps'] != 4608 or
            value['revision' if state else 'feature_revision'] != F.A.R.C.I.F.REVISION):
        raise ValueError('Only the fixed complete A37 last004608 model is allowed')
    if not state and (value['completed'] is not True or value['preflight_job'] != F.PREFLIGHT_JOB or
            value['trainable_parameters'] != 612324 or value['initialization'] != CONFIG['initialization'] or
            value['final_checkpoint'] != str((F.OUT/'step_004608.pt').relative_to(ROOT)) or
            value['reference_injected_training_candidates'] != 0 or
            any(value[k] is not False for k in ('supervised_checkpoint_read','generator_forward',
                'generator_updated','new_candidate_search','new_backbone_forward','new_R4_labelling',
                'DEV_used_for_selection','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual complete unchanged natural-only A37 fit required')


def gate(counts):
    checks = {}
    for species, values in counts.items():
        old, new = values['A31_A3'], values['A37']
        checks[species+'_natural_selected_correct_strictly_increased'] = new['selected_correct'] > old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = (
            new['selected_known_nonmatch'] <= old['selected_known_nonmatch'])
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A37']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def paired_readout(saved, previous, labelled):
    if not len(saved) == len(previous) == len(labelled):
        raise ValueError('Complete paired original windows required')
    totals = {s: {a: Counter() for a in ('A31_A3','A37')} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    rows = []
    for new, old, lab in zip(saved, previous, labelled):
        if (new['window'] != old['window'] or new['window'] != lab['window'] or
                new['reference_used_in_scoring'] or new['chains'] != old['variants']['A31_A3']['chains']):
            raise ValueError('Unchanged ordered reference-free original natural geometry required')
        labels = lab['labels']['A31_A3']
        if any(y not in (-1,0,1) for y in labels):
            raise ValueError('Original protected labels required')
        counts = dict(A31_A3=R.selected_counts(old['variants']['A31_A3'], labels),
                      A37=R.selected_counts(new, labels))
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
    draws, free, previous, index, old = B.population()
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    baseline = [dict(window=r['window'], reference_used_in_scoring=False,
                     **r['variants']['A31_A3']) for r in previous]
    rows, counts, references = paired_readout(baseline, previous, labelled)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A37'] != counts[s]['A31_A3'] for s in counts)):
        raise ValueError('Original metadata pairing and denominators differ')
    R.C.R.write_summary(CHECKS, dict(completed=True, experiment='M28-A37-CONTEXT-CHAIN-READOUT-CHECKS',
        protocol=PROTOCOL, configuration=CONFIG, windows=1491, original_draws=len(draws),
        natural_candidates=sum(len(v['chains']) for v in free), feature_index_records=len(index),
        original_identity_labels_DP_stronger_baseline_and_denominators_reconciled=True,
        new_fitted_checkpoint_opened=False, feature_tensors_opened=False, model_or_backbone_forward=False,
        optimizer_steps_this_run=0, new_R4_labelling=False, new_candidate_search=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Original saved metadata plus synthetic checks only, not new-model performance.'))


@torch.inference_mode()
def main(checks_job):
    F.A.require_jobs([FIT_JOB, str(checks_job)])
    ready = R.completed(CHECKS, checks_job)
    if (ready['protocol'] != PROTOCOL or ready['configuration'] != CONFIG or
            ready['original_identity_labels_DP_stronger_baseline_and_denominators_reconciled'] is not True):
        raise ValueError('Actual matching readout readiness required')
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
    rows, counts, references = paired_readout(saved, previous, labelled)
    if (len(rows) != 1491 or sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Original full denominators and stronger baseline differ')
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    checks = gate(counts)
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A37-CONTEXT-CHAIN-TRAIN-READOUT',
        protocol=PROTOCOL, fit_protocol=F.PROTOCOL, configuration=CONFIG, source_fit_job=FIT_JOB,
        source_checks_job=str(checks_job), source_baseline_jobs=dict(readout='13398521', CPU='13398524'),
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
    parser.add_argument('--mode', choices=('checks','readout'), required=True)
    parser.add_argument('--checks-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'checks':
        readiness()
    else:
        if not args.checks_job:
            parser.error('--checks-job is required for actual readout')
        main(args.checks_job)
