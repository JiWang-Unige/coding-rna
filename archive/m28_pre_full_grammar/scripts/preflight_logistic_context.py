#!/usr/bin/env python3
"""A38 zero-update pointwise-null objective in the unchanged joint model."""
import math
import time
import torch
import preflight_context_chain as A
import readout_context_chain_train as D
from src.m28.core import chain_loss

ROOT = A.ROOT
R, B = A.R, A.B
OUT = ROOT/'outputs/M28-A38-LOGISTIC-CONTEXT-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a38_logistic_context_protocol.md'
CONFIG = {**A.CONFIG, 'loss': 'mean_known_natural_candidate_BCE_against_window_null',
          'chain_candidate_BCE': True,
          'loss_unknown': 'exclude_-1_from_BCE_keep_in_natural_competition'}
new_model, population, require_jobs = A.new_model, A.population, A.require_jobs


def closed_a37():
    require_jobs(['13410061','13410739','13410742','13412326'])
    fitted = A.B.R.completed(D.F.OUT, '13410061'); D.validate_fit(fitted)
    report = A.B.R.completed(D.OUT, '13410739')
    checked = A.B.R.completed(D.CHECK, '13410742')
    transition = A.B.R.completed(ROOT/'outputs/M28-A37-CONTEXT-CHAIN-ATTRITION', '13412326')
    proofs = ('actual_original_draw_weights_counts_and_4608_updates_verified',
              'actual_three_model_checkpoint_tensors_and_branch_optimizer_exposure_verified',
              'actual_saved_natural_identity_and_FP32_score_components_verified',
              'actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed')
    if (report['configuration'] != A.CONFIG or checked['configuration'] != A.CONFIG or
            checked['source_jobs'] != dict(fit='13410061',readout='13410739') or
            any(checked[k] is not True for k in proofs) or
            checked['selected_counts_per_species'] != report['selected_counts_per_species'] or
            checked['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            transition['selected_counts_per_species'] != report['selected_counts_per_species'] or
            transition['windows'] != 1491 or transition['natural_candidates'] != 16321 or
            report['TRAIN_preparation_supported'] is not False or
            checked['TRAIN_preparation_supported'] is not False or
            D.gate(report['selected_counts_per_species']) != report['predeclared_TRAIN_gate']):
        raise ValueError('Closed numerically verified A37 required; this is not its continuation')


def main():
    closed_a37()
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
        loss = chain_loss(score, labels); weighted = loss*weights[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual mean-known BCE')
        score['logits'].retain_grad(); weighted.backward()
        gain_grad = score['logits'].grad
        known = sum(y >= 0 for y in labels)
        expected = torch.zeros_like(score['logits'])
        for i, y in enumerate(labels):
            if y >= 0:
                expected[i] = (torch.sigmoid(score['logits'][i].detach())-y)*weights[ident]/known
            elif float(gain_grad[i]) != 0.:
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
                raise ValueError('Missing/nonfinite actual BCE branch gradient: '+name)
            norms[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if norms[name] <= 0:
                raise ValueError('Zero applicable branch gradient: '+name)
        torch.cuda.synchronize()
        receipts.append(dict(**witness, weighted_loss=float(weighted.detach()), gradient_norms=norms,
            known_junction_branch_applicable=applicable, unknown_direct_gain_gradient_zero=True,
            actual_mean_BCE_gain_gradient_verified=True,
            load_forward_backward_seconds=time.perf_counter()-started,
            known_candidates=known, known_positives=labels.count(1),
            known_nonmatches=labels.count(0), unknown_excluded=labels.count(-1)))
        print(receipts[-1], flush=True)
        del cached, features, dna, score, loss, weighted, gain_grad, expected
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            not all(any(r['species'] == s and r['known_junction_branch_applicable'] for r in receipts)
                    for s in A.R.C.R.ALLOW)):
        raise ValueError('Unchanged parameters and known junction coverage in both species required')
    A.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A38-LOGISTIC-CONTEXT-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, windows=1491, original_draws=1536, natural_candidates=16321,
        natural_counts=A.R.totals(records, records),
        replayed_original_draw_counts=A.R.totals(records, [d['window_id'] for d in draws]),
        eligible_reference_pairs=references, stronger_baseline_counts=baseline['selected_counts_per_species'],
        original_natural_identity_labels_and_stronger_baseline_reconciled=True,
        all_natural_gold_sets_compatible=True, fresh_weights_unchanged=True,
        both_species_known_junction_gradient_coverage=True, actual_mean_BCE_gain_gradient_verified=True,
        initialization_identical_constructor_and_seed_to_A37=True, original_TRAIN_catalog_metadata_read=True,
        supervised_checkpoint_read=False, generator_forward=False, optimizer_steps_this_run=0,
        checkpoints_written=False, new_candidate_search=False, new_backbone_forward=False,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        reference_injected_training_candidates=0, records=receipts, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update pointwise-null objective readiness; not trained performance or scientific PASS.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
