#!/usr/bin/env python3
"""Scalar mechanism evidence only: no research data, model or optimizer."""
import argparse
import itertools
import json
import math
from pathlib import Path

import torch
from src.m28.core import chain_loss
from src.m28.chain_set import chain_set_loss
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss


def compatible(subset, chains):
    ordered = sorted(subset, key=lambda i: chains[i][0][0])
    return all(chains[a][-1][1] <= chains[b][0][0]
               for a, b in zip(ordered, ordered[1:]))


def enumerate_sets(chains):
    return [tuple(i for i in range(len(chains)) if mask & (1 << i))
            for mask in range(1 << len(chains))
            if compatible(tuple(i for i in range(len(chains)) if mask & (1 << i)), chains)]


def enumerate_partition(values, subsets):
    weights = [math.exp(sum(values[i] for i in subset)) for subset in subsets]
    partition = sum(weights)
    occupancy = [sum(w for w, s in zip(weights, subsets) if i in s)/partition
                 for i in range(len(values))]
    return math.log(partition), occupancy


def close(actual, expected):
    torch.testing.assert_close(actual, actual.new_tensor(expected), rtol=1e-11, atol=1e-11)


def verify_case(values, chains, gold):
    # The high-scoring unknown overlaps every known chain but is not supervised.
    all_chains = chains + [[[0, 12]]]
    labels = [int(i in gold) for i in range(len(chains))] + [-1]
    gain = torch.tensor(list(values)+[7.], dtype=torch.float64, requires_grad=True)
    scores = dict(logits=gain, null=gain.new_tensor(.2, requires_grad=True))
    bce = chain_loss(scores, labels)
    structured = chain_set_loss(scores, all_chains, labels)['loss']
    composite = null_anchored_chain_set_loss(scores, all_chains, labels)['loss']
    logz, occupancy = enumerate_partition(values, enumerate_sets(chains))
    sigmoid = [1/(1+math.exp(-v)) for v in values]
    log_all = sum(math.log1p(math.exp(v)) for v in values)
    n = len(chains)
    gb, gs, gc = [torch.autograd.grad(loss, gain, retain_graph=True)[0]
                  for loss in (bce, structured, composite)]
    close(structured, float(bce.detach())+(logz-log_all)/n)
    close(composite, float(bce.detach())+(logz-log_all)/(2*n))
    close(gb, [(p-y)/n for p, y in zip(sigmoid, labels)]+[0.])
    close(gs, [(q-y)/n for q, y in zip(occupancy, labels)]+[0.])
    close(gc, [((p+q)/2-y)/n for p, q, y in zip(sigmoid, occupancy, labels)]+[0.])
    close(gc-gb, [(q-p)/(2*n) for p, q in zip(sigmoid, occupancy)]+[0.])
    if any(q > p+1e-11 for q, p in zip(occupancy, sigmoid)):
        raise AssertionError('Structured occupancy exceeds independent Bernoulli occupancy')
    if bool((gc > gb+1e-11).any()) or bool((gc < gs-1e-11).any()):
        raise AssertionError('Composite independent-gain derivative ordering differs')
    return dict(gain=list(values), labels=labels[:-1], compatible_sets=len(enumerate_sets(chains)),
                known=len(chains), sigmoid=sigmoid, structured_occupancy=occupancy,
                log_compatibility_probability=logz-log_all,
                losses=dict(BCE=float(bce.detach()), set_NLL=float(structured.detach()),
                            composite=float(composite.detach())),
                gradients=dict(BCE=gb.tolist(), set_NLL=gs.tolist(), composite=gc.tolist()))


def shared_null_example():
    chains = [[[0, 10]], [[1, 9]], [[2, 8]], [[0, 12]]]
    labels = [1, 0, 0, -1]
    raw = torch.tensor([8.4, 2.4, -1.6, 7.4], dtype=torch.float64, requires_grad=True)
    null = torch.tensor(.4, dtype=torch.float64, requires_grad=True)
    scores = dict(logits=raw-null, null=null)
    losses = (chain_loss(scores, labels), chain_set_loss(scores, chains, labels)['loss'],
              null_anchored_chain_set_loss(scores, chains, labels)['loss'])
    gradients = [torch.autograd.grad(loss, (raw, null), retain_graph=True) for loss in losses]
    gb, gs, gc = [g[0] for g in gradients]
    nb, ns, nc = [g[1] for g in gradients]
    close(nc-nb, -float((gc-gb).sum()))
    close(nc, (float(nb)+float(ns))/2)
    if not float(nc-nb) > 0 or any(float(g[-1]) != 0. for g in (gb, gs, gc)):
        raise AssertionError('Unknown direct raw-energy exclusion or shared-null shift differs')
    return dict(gain=scores['logits'].detach().tolist(), labels=labels,
                null_gradients=dict(BCE=float(nb), set_NLL=float(ns), composite=float(nc)),
                composite_minus_BCE_null_gradient=float(nc-nb), unknown_raw_direct_gradient=0.,
                unknown_gain_derivative_wrt_null=-1.,
                limitation='Independent raw-energy/null scalar coordinates; not fitted neural parameter dynamics')


def main(out):
    layouts = dict(disjoint=[[[0, 3]], [[3, 6]], [[6, 9]]],
                   clique=[[[0, 9]], [[1, 8]], [[2, 7]]],
                   path=[[[0, 4]], [[3, 7]], [[6, 10]]],
                   pair_and_isolated=[[[0, 5]], [[1, 4]], [[6, 10]]])
    counts = {}
    for name, chains in layouts.items():
        subsets = enumerate_sets(chains)
        count = 0
        for values in itertools.product((-3., 0., 3.), repeat=3):
            for gold in subsets:
                verify_case(values, chains, gold)
                count += 1
        counts[name] = count
    if sum(counts.values()) != 621:
        raise AssertionError('The fixed synthetic design differs')
    examples = dict(strong_positive_conflict=verify_case((8., 2.), [[[0, 10]], [[1, 9]]], (0,)),
                    losing_positive_conflict=verify_case((.5, 1.), [[[0, 10]], [[1, 9]]], (0,)),
                    compatible_multiple_genes=verify_case((.5, 1.), [[[0, 5]], [[5, 10]]], (0, 1)))
    report = dict(experiment='M28-A39-LOSS-GRADIENT-MECHANISM', completed=True,
                  scope='deterministic_synthetic_scalar_mechanism_not_research_model_performance',
                  fixed_gain_grid=[-3., 0., 3.], layouts=layouts, cases_per_layout=counts,
                  total_labelled_grid_cases=621, examples=examples,
                  shared_null_example=shared_null_example(), arithmetic_dtype='float64',
                  original_loss_implementations_unchanged=True,
                  conditional_Bernoulli_identity_verified=True,
                  exact_enumeration_derivative_and_occupancy_bound_verified=True,
                  label_independent_residual_derivative_verified=True,
                  shared_null_and_direct_unknown_distinction_verified=True,
                  research_data_used=False, model_constructed=False, checkpoint_used=False,
                  optimizer_steps=0, backbone_forward=False, DEV_used=False,
                  test_setaria_used=False, A33_outputs_opened=False,
                  limitation='Universal inequality follows from the analytic partition proof, not this finite grid. '
                  'Direct scalar derivatives do not predict shared neural weights, AdamW, selected chains or F1.')
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True, cases_per_layout=counts,
                         total_labelled_grid_cases=621, optimizer_steps=0), allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.output_dir)
