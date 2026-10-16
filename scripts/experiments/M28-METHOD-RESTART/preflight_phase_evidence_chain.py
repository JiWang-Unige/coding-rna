#!/usr/bin/env python3
"""A26 zero-update engineering readiness on the original eight TRAIN cases."""
import math
import time
import torch
import run_junction_chain_set as J
from src.m28.phase_evidence_chain import PhaseEvidenceChainHead, phase_evidence

R, ROOT = J.R, J.ROOT
OUT = ROOT/'outputs/M28-A26-PHASE-EVIDENCE-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a26_phase_evidence_protocol.md'
CONFIG = {**J.CONFIG, 'quality_head': 'PhaseEvidenceChainHead_A21_plus_frozen_state_odds',
    'quality_trainable_parameters': 248096, 'new_parameters': 28,
    'phase_evidence': 'state1_to14_mean_logit_difference_vs0_then_log1p_role_counts',
    'phase_projection': 'linear28_to1_no_bias_zero_initialization',
    'phase_local_predictor': 'frozen_A3_segmentation_same_last004608',
    'main_evaluation': 'last_head_step_004608_only_natural_TRAIN_A6_gate_A21_pair'}
GROUPS = J.GROUPS+('path_evidence',)


def new_head(device):
    torch.manual_seed(CONFIG['quality_seed'])
    head = PhaseEvidenceChainHead(CONFIG['hidden']).to(device)
    if sum(p.numel() for p in head.parameters()) != CONFIG['quality_trainable_parameters']:
        raise ValueError('Exact A26 head structure required')
    return head


@torch.no_grad()
def evidence_inputs(generator, hidden, window, record):
    chains = [row['chain'] for row in record['candidates']]
    evidence = phase_evidence(generator.segmentation(hidden), chains, window['valid_bases'])
    if evidence.requires_grad or not bool(torch.isfinite(evidence).all()):
        raise ValueError('Expected frozen finite candidate-state evidence')
    return evidence


def score_loss(head, hidden, evidence, window, record):
    rows = record['candidates']; chains = [row['chain'] for row in rows]
    scores = head(hidden, chains, evidence, window['valid_bases'])
    result = R.chain_set_loss(scores, chains, [row['label'] for row in rows])
    if any(not bool(torch.isfinite(value).all()) for value in scores.values()) or not bool(torch.isfinite(result['loss'])):
        raise ValueError('Nonfinite A26 actual score/structured loss')
    return scores, result


def zero_update_witnesses(head, generator, index, windows, records, hidden_by_id=None, evidence_by_id=None):
    initial = {k: v.detach().clone() for k, v in head.state_dict().items()}
    receipts = []
    for witness in R.P.select_gradient_windows(records):
        ident = witness['window_id']; window = windows[ident]; local = records[ident]
        chains = [row['chain'] for row in local['candidates']]
        junction_applicable = any(row['label'] >= 0 and len(row['chain']) > 1 for row in local['candidates'])
        live = R.P.inputs(generator, index, window).detach()
        evidence = evidence_inputs(generator, live, window, local)
        values = ((live, evidence),) if hidden_by_id is None else ((live, evidence), (hidden_by_id[ident], evidence_by_id[ident]))
        if hidden_by_id is not None:
            torch.testing.assert_close(values[0], values[1], rtol=1e-6, atol=1e-6)
        prior = None
        for hidden, frozen_evidence in values:
            head.zero_grad(set_to_none=True)
            scores, result = score_loss(head, hidden, frozen_evidence, window, local)
            scores['logits'].retain_grad(); result['loss'].backward()
            if any(float(scores['logits'].grad[i]) != 0. for i, row in enumerate(local['candidates']) if row['label'] == -1):
                raise ValueError('Unknown candidate received direct gain supervision')
            groups = {}
            for name in GROUPS:
                params = list(getattr(head, name).parameters())
                if name in ('junction', 'junction_additive') and not junction_applicable:
                    if any(p.grad is not None and bool((p.grad != 0).any()) for p in params):
                        raise ValueError('Unsupported junction branch received gradient')
                    groups[name] = 0.
                    continue
                if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                    raise ValueError('Missing/nonfinite actual gradient: '+name)
                groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
                if groups[name] <= 0:
                    raise ValueError('Inactive actual gradient branch: '+name)
            norm = torch.nn.utils.clip_grad_norm_(head.parameters(), CONFIG['clip_norm'], error_if_nonfinite=True)
            snapshot = dict(scores={k: v.detach().clone() for k, v in scores.items()},
                loss=result['loss'].detach().clone(), norm=norm.detach().clone(),
                gradient_present=torch.tensor([p.grad is not None for p in head.parameters()]),
                gradients={k: p.grad.detach().clone() if p.grad is not None else torch.zeros_like(p)
                           for k, p in head.named_parameters()})
            if prior is not None:
                torch.testing.assert_close(snapshot, prior, rtol=1e-6, atol=1e-6)
            prior = snapshot
        if any(p.grad is not None for p in generator.parameters()):
            raise ValueError('Frozen generator or local predictor received gradient')
        receipts.append(dict(**witness, gradient_norms=groups, loss=float(result['loss'].detach()),
            candidate_rows=len(chains), evidence_shape=list(evidence.shape), evidence_finite=True,
            junction_branch_applicable=junction_applicable, unknown_direct_gain_gradient_zero=True,
            reference_used_in_evidence=False, live_cached_equivalence=hidden_by_id is not None))
    if any(not torch.equal(initial[k], v) for k, v in head.state_dict().items()):
        raise ValueError('Zero-update witness check changed head weights')
    if any(not any(row['species'] == species and row['gradient_norms']['path_evidence'] > 0 for row in receipts)
           for species in R.C.R.ALLOW):
        raise ValueError('Both species require actual new evidence projection gradient')
    head.zero_grad(set_to_none=True)
    return receipts


def main():
    began = time.perf_counter()
    draws, index, windows, records, _ = R.population()
    if len(draws) != 1536 or len(records) != 1491:
        raise ValueError('Original TRAIN population required')
    generator = R.C.load_generator(); head = new_head('cuda').train()
    original = J.new_head('cuda').eval().requires_grad_(False)
    for key, value in original.state_dict().items():
        if not torch.equal(value, head.state_dict()[key]):
            raise ValueError('Fresh A21 common parameter mismatch: '+key)
    comparisons = []
    with torch.no_grad():
        for witness in R.P.select_gradient_windows(records):
            ident = witness['window_id']; window = windows[ident]; local = records[ident]
            hidden = R.P.inputs(generator, index, window).detach()
            evidence = evidence_inputs(generator, hidden, window, local)
            current, loss = score_loss(head, hidden, evidence, window, local)
            prior, old_loss = J.score_loss(original, hidden, window, local)
            for key, value in prior.items():
                torch.testing.assert_close(value, current[key], rtol=0, atol=0)
            torch.testing.assert_close(loss['loss'], old_loss['loss'], rtol=0, atol=0)
            comparisons.append(dict(**witness, zero_prior_A21_scores_and_loss_match=True))
    del original
    receipts = zero_update_witnesses(head, generator, index, windows, records)
    R.C.R.write_summary(OUT, dict(completed=True, arm='A26', configuration=CONFIG, protocol=PROTOCOL,
        original_draws=1536, windows=1491, natural_counts=R.totals(records, records),
        replayed_draw_counts=R.totals(records, [draw['window_id'] for draw in draws]),
        quality_trainable_parameters=248096, exact_A21_common_fresh_parameters=True,
        zero_prior_actual_scores_and_loss_verified=True, comparisons=comparisons,
        gradient_witnesses=receipts, head_weights_unchanged=True, generator_frozen=True,
        optimizer_steps_this_run=0, new_candidate_search=False, reference_used_in_evidence=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='Zero-update engineering readiness only; no TRAIN selection or independent accuracy result.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
