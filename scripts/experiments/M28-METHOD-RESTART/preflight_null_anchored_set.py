#!/usr/bin/env python3
"""A39 zero-update absolute-null plus compatible-set objective in the joint model."""
import math
import time
import torch
import preflight_context_chain as A
import readout_logistic_context_train as L
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss

ROOT = A.ROOT
R, B = A.R, A.B
OUT = ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-PREFLIGHT-R1'
PROTOCOL = 'reports/M28-METHOD-RESTART/a39_null_anchored_set_protocol.md'
CONFIG = {**A.CONFIG, 'loss': 'equal_mean_known_BCE_and_compatible_set_NLL',
          'chain_candidate_BCE': True, 'compatible_known_set_NLL': True,
          'fixed_objective_weights': [.5, .5],
          'loss_unknown': 'exclude_-1_from_BCE_partition_and_gold_keep_in_natural_competition'}
new_model, population, require_jobs = A.new_model, A.population, A.require_jobs


SCALAR_JOB = '13415040'
CLOSED_A38_JOBS = ['13412630', '13412701', '13412702', '13414584', '13414782']


def closed_a38():
    require_jobs(CLOSED_A38_JOBS+[SCALAR_JOB])
    scalar = ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-CHECKS'
    if (scalar.joinpath('STATUS').read_text().strip() != 'COMPLETED' or
            scalar.joinpath('JOBID').read_text().strip() != SCALAR_JOB):
        raise ValueError('Completed actual synthetic A39 objective checks required')
    fitted = A.B.R.completed(L.F.OUT, '13412630'); L.validate_fit(fitted)
    report = A.B.R.completed(L.OUT, '13412701')
    checked = A.B.R.completed(L.CHECK, '13412702')
    attrition = A.B.R.completed(ROOT/'outputs/M28-A38-LOGISTIC-CONTEXT-ATTRITION', '13414584')
    competition = A.B.R.completed(ROOT/'outputs/M28-A38-LOGISTIC-CONTEXT-COMPETITION', '13414782')
    proofs = ('actual_original_draw_weights_counts_and_4608_updates_verified',
              'actual_three_model_checkpoint_tensors_and_branch_optimizer_exposure_verified',
              'actual_saved_natural_identity_and_FP32_score_components_verified',
              'actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed')
    if (report['configuration'] != L.CONFIG or checked['configuration'] != L.CONFIG or
            checked['source_jobs'] != dict(fit='13412630', readout='13412701') or
            any(checked[k] is not True for k in proofs) or
            checked['selected_counts_per_species'] != report['selected_counts_per_species'] or
            checked['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            report['TRAIN_preparation_supported'] is not False or
            checked['TRAIN_preparation_supported'] is not False or
            L.gate(report['selected_counts_per_species']) != report['predeclared_TRAIN_gate'] or
            competition['all_individual_and_joint_actual_selected_exchanges_verified'] is not True or
            competition['source_jobs'] != dict(fit='13412630', readout='13412701', CPU='13412702',
                attrition='13414584', baseline_readout='13398521')):
        raise ValueError('Closed numerically verified A38 required; this is not its continuation')
    for value in (attrition, competition):
        if (value['configuration'] != L.CONFIG or value['windows'] != 1491 or
                value['natural_candidates'] != 16321 or
                value['scientific_status'] != 'A38_TRAIN_FAIL_closed_no_rescue' or
                {s: {arm: counts[arm] for arm in ('A31_A3', 'A38')}
                 for s, counts in value['selected_counts_per_species'].items()} !=
                {s: {arm: counts[arm] for arm in ('A31_A3', 'A38')}
                 for s, counts in report['selected_counts_per_species'].items()} or
                value['eligible_reference_pairs'] != report['eligible_reference_pairs']):
            raise ValueError('Actual original closed A38 diagnosis required')


def main():
    closed_a38()
    draws, index, windows, records, source, baseline, references = population()
    model = new_model('cuda').train()
    if sum(p.numel() for p in model.parameters()) != CONFIG['trainable_parameters']:
        raise ValueError('Unchanged fresh A37 architecture required')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    weights = {d['window_id']: d['importance_weight_uniform_within_stratum'] for d in draws}
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); receipts = []
    for witness in A.R.P.select_gradient_windows(records):
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        chains, labels = [r['chain'] for r in local], [r['label'] for r in local]
        applicable = any(y >= 0 and len(c) > 1 for c, y in zip(chains, labels))
        started = time.perf_counter()
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        features, dna = A.R.C.R.Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Original encoder features and DNA must remain detached')
        model.zero_grad(set_to_none=True)
        score = model(features, dna, [w['valid_bases']], chains)
        if any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
            raise ValueError('Finite unchanged FP32 complete-chain scores required')
        result = null_anchored_chain_set_loss(score, chains, labels)
        loss = result['loss']; weighted = loss*weights[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual fixed composite objective')
        BCE_grad = torch.autograd.grad(result['mean_known_BCE'], score['logits'], retain_graph=True)[0]
        set_grad = torch.autograd.grad(result['mean_known_set_NLL'], score['logits'], retain_graph=True)[0]
        expected = (BCE_grad+set_grad)*weights[ident]/2
        score['logits'].retain_grad(); weighted.backward()
        gain_grad = score['logits'].grad
        known = result['known_candidates']
        for i, y in enumerate(labels):
            if y == -1 and float(gain_grad[i]) != 0.:
                raise ValueError('Unknown must have exactly zero direct gain gradient')
        torch.testing.assert_close(gain_grad, expected, rtol=1e-6, atol=1e-7)
        modules = [('shared.'+k, getattr(model.shared, k)) for k in ('project','context','dna','fuse')]
        modules += [('head.'+k, getattr(model.head, k)) for k in A.HEAD_GROUPS]
        norms = {}
        for name, module in modules:
            params = list(module.parameters())
            if name in ('head.junction','head.junction_additive') and not applicable:
                if any(p.grad is not None and not bool((p.grad == 0).all()) for p in params):
                    raise ValueError('No known junction must give no junction gradient: '+name)
                norms[name] = 0.; continue
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual composite branch gradient: '+name)
            norms[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if norms[name] <= 0:
                raise ValueError('Zero applicable branch gradient: '+name)
        torch.cuda.synchronize()
        receipts.append(dict(**witness, weighted_loss=float(weighted.detach()),
            mean_known_BCE=float(result['mean_known_BCE'].detach()),
            mean_known_set_NLL=float(result['mean_known_set_NLL'].detach()), gradient_norms=norms,
            known_junction_branch_applicable=applicable, unknown_direct_gain_gradient_zero=True,
            actual_fixed_composite_gain_gradient_verified=True,
            load_forward_backward_seconds=time.perf_counter()-started,
            known_candidates=known, known_positives=labels.count(1),
            known_nonmatches=labels.count(0), unknown_excluded=labels.count(-1)))
        print(receipts[-1], flush=True)
        del cached, features, dna, score, result, loss, weighted, gain_grad, expected, BCE_grad, set_grad
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            not all(any(r['species'] == s and r['known_junction_branch_applicable'] for r in receipts)
                    for s in A.R.C.R.ALLOW)):
        raise ValueError('Unchanged parameters and known junction coverage in both species required')
    A.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A39-NULL-ANCHORED-SET-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, source_scalar_checks_job=SCALAR_JOB,
        source_closed_A38_jobs=CLOSED_A38_JOBS, windows=1491, original_draws=1536, natural_candidates=16321,
        natural_counts=A.R.totals(records, records),
        replayed_original_draw_counts=A.R.totals(records, [d['window_id'] for d in draws]),
        eligible_reference_pairs=references, stronger_baseline_counts=baseline['selected_counts_per_species'],
        original_natural_identity_labels_and_stronger_baseline_reconciled=True,
        all_natural_gold_sets_compatible=True, fresh_weights_unchanged=True,
        both_species_known_junction_gradient_coverage=True, actual_fixed_composite_gain_gradient_verified=True,
        initialization_identical_constructor_and_seed_to_A37=True, original_TRAIN_catalog_metadata_read=True,
        supervised_checkpoint_read=False, generator_forward=False, optimizer_steps_this_run=0,
        checkpoints_written=False, new_candidate_search=False, new_backbone_forward=False,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        reference_injected_training_candidates=0, records=receipts, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update fixed absolute-null plus compatible-set readiness; not trained performance or scientific PASS.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
