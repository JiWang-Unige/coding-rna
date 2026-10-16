"""A training-only one-to-one query assignment; never used by DNA inference.

All positive chains must fit. Unmatched labels and loss coefficients are
explicit caller decisions, not inferred from missing/partial annotations.
Rectangular linear assignment is a discrete supervision step, not a novel
decoder or a claim of a globally differentiable annotation system.
"""
import math
import torch
from torch.nn import functional as F
from scipy.optimize import linear_sum_assignment


@torch.no_grad()
def teacher_chain_components(window, chain):
    """Evaluate one reference chain for every query, only for assignment costs.

    Uses the raw teacher's geometry and event means, not the DNA decoder's
    minimum-intron rules. No states are retained for gradient recomputation.
    """
    window.validate_chain(chain)
    states, previous_end = window.initial, -1
    q = window.model.num_queries
    boundaries, endings = [], []
    for j, (a, b) in enumerate(chain):
        stage = 'start' if j == 0 else 'next_start'
        lower = 0 if j == 0 else previous_end+1
        boundaries.append(F.cross_entropy(
            window.pointer_logits(states, stage, lower),
            torch.full((q,), a, device=states.device, dtype=torch.long),
            reduction='none'))
        states = window.advance(states, stage, a)
        boundaries.append(F.cross_entropy(
            window.pointer_logits(states, 'end', a),
            torch.full((q,), b, device=states.device, dtype=torch.long),
            reduction='none'))
        states = window.advance(states, 'end', b)
        logits = window.model.terminate(states).squeeze(-1)
        endings.append(F.binary_cross_entropy_with_logits(
            logits, logits.new_full((q,), float(j == len(chain)-1)),
            reduction='none'))
        previous_end = b
    return {'boundary_ce': torch.stack(boundaries).mean(dim=0),
            'termination_bce': torch.stack(endings).mean(dim=0)}


def matched_query_components(model, features, positive_chains, coefficients,
                             unmatched_label, valid_bases=None):
    """Match all supplied positive chains, then recompute a gradient-bearing loss.

    coefficients order: existence BCE, mean boundary CE, mean termination BCE.
    unmatched_label is 0 (explicit reference-relative negative policy) or -1
    (unknown/protected); no default labels/coefficients or genome files here.
    Caller still owns R4 eligibility/alternative/partial masks and capacity.
    """
    coefficients = tuple(float(c) for c in coefficients)
    if len(coefficients) != 3 or any(not math.isfinite(c) or c < 0 for c in coefficients):
        raise ValueError('Three finite nonnegative loss coefficients are required')
    if unmatched_label not in (-1, 0):
        raise ValueError('Unmatched queries require explicit unknown or negative labels')
    chains = tuple(tuple(tuple(exon) for exon in chain) for chain in positive_chains)
    n, q = len(chains), model.num_queries
    if n > q:
        raise ValueError('Query capacity cannot drop supplied positive chains')
    window = model.for_window(features, valid_bases)
    with torch.no_grad():
        costs = window.query_logits.new_empty((0, q))
        gene_to_query = ()
        if n:
            # The all-negative existence term is a constant when unmatched
            # labels are zero. Assignment cost uses positive-minus-negative,
            # not positive BCE alone; normalizers match teacher_components.
            known = q if unmatched_label == 0 else n
            delta = F.softplus(-window.query_logits)
            if unmatched_label == 0:
                delta = delta - F.softplus(window.query_logits)
            rows = []
            for chain in chains:
                pair = teacher_chain_components(window, chain)
                rows.append(coefficients[0]*delta/known +
                            coefficients[1]*pair['boundary_ce']/n +
                            coefficients[2]*pair['termination_bce']/n)
            costs = torch.stack(rows)
            if not torch.isfinite(costs).all():
                raise ValueError('Nonfinite query assignment cost')
            gene_indices, query_indices = linear_sum_assignment(costs.cpu().numpy())
            assignment = [None]*n
            for gene, query in zip(gene_indices, query_indices):
                assignment[int(gene)] = int(query)
            gene_to_query = tuple(assignment)

    labels = [unmatched_label]*q
    assigned_chains = {}
    for gene, query in enumerate(gene_to_query):
        labels[query] = 1
        assigned_chains[query] = chains[gene]
    # No no-grad hidden state or cost graph is used for the trainable readout.
    components = window.teacher_components(labels, assigned_chains)
    names = ('existence_bce', 'boundary_ce', 'termination_bce')
    total = sum(c*components[name] for c, name in zip(coefficients, names))
    return {'total': total, 'components': components, 'query_labels': tuple(labels),
            'positive_chains': assigned_chains, 'gene_to_query': gene_to_query,
            'assignment_costs': costs}
