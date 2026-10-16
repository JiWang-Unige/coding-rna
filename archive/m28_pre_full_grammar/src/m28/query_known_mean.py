"""A3 known-query mean existence loss; A2 architecture and inference are unchanged."""
import math
import torch
from torch.nn import functional as F
from scipy.optimize import linear_sum_assignment
from .labels import local_loss
from .query_model import unmatched_query_policy
from .query_training import teacher_chain_components


def known_query_components(window, positive_chains, coefficients, unmatched_label):
    """Assign every positive using the exact known-query-mean window objective.

    Unknown -1 logits receive no direct BCE. Boundary and EOS means are exactly
    the existing teacher components. No reference or labels enter free inference.
    """
    coefficients = tuple(float(c) for c in coefficients)
    if len(coefficients) != 3 or any(not math.isfinite(c) or c < 0 for c in coefficients):
        raise ValueError('Three finite nonnegative query coefficients are required')
    if unmatched_label not in (-1, 0):
        raise ValueError('Explicit unknown or known-negative unmatched label required')
    chains = tuple(tuple(tuple(exon) for exon in chain) for chain in positive_chains)
    n, q = len(chains), window.model.num_queries
    if n > q:
        raise ValueError('Query capacity cannot drop supplied positive chains')
    with torch.no_grad():
        costs = window.query_logits.new_empty((0, q))
        gene_to_query = ()
        if n:
            known = q if unmatched_label == 0 else n
            delta = F.softplus(-window.query_logits)
            if unmatched_label == 0:
                # All-negative BCE is assignment-independent. Replacing each
                # matched label with 1 changes BCE by positive minus negative.
                delta = delta - F.softplus(window.query_logits)
            rows = []
            for chain in chains:
                pair = teacher_chain_components(window, chain)
                rows.append(coefficients[0]*delta/known +
                            coefficients[1]*pair['boundary_ce']/n +
                            coefficients[2]*pair['termination_bce']/n)
            costs = torch.stack(rows)
            if not torch.isfinite(costs).all():
                raise ValueError('Nonfinite A3 query assignment cost')
            genes, queries = linear_sum_assignment(costs.cpu().numpy())
            assignment = [None]*n
            for gene, query in zip(genes, queries):
                assignment[int(gene)] = int(query)
            gene_to_query = tuple(assignment)
    labels = [unmatched_label]*q
    assigned = {}
    for gene, query in enumerate(gene_to_query):
        labels[query] = 1
        assigned[query] = chains[gene]
    # Existing teacher averages BCE over all known labels, not over classes.
    # Recompute with gradients; never reuse the discrete matching graph.
    components = window.teacher_components(labels, assigned)
    names = ('existence_bce', 'boundary_ce', 'termination_bce')
    return {'total': sum(c*components[name] for c, name in zip(coefficients, names)),
            'components': components, 'query_labels': tuple(labels),
            'positive_chains': assigned, 'gene_to_query': gene_to_query,
            'assignment_costs': costs}


def joint_known_query_losses(model, outputs, targets, valid_bases, window_weight, coefficients):
    """Use the unchanged local labels and unknown policy with the A3 objective."""
    coefficients = tuple(float(c) for c in coefficients)
    weight = float(window_weight)
    if len(coefficients) != 4 or any(not math.isfinite(c) or c < 0 for c in coefficients):
        raise ValueError('Four finite nonnegative joint coefficients are required')
    if not math.isfinite(weight) or weight < 0:
        raise ValueError('Finite nonnegative window importance weight required')
    if outputs['features'].shape[0] != 1:
        raise ValueError('Expected a single-window query loss')
    window = model.query_window(outputs['features'][0], valid_bases)
    policy = unmatched_query_policy(targets, valid_bases)
    query = known_query_components(window, targets['positive_chains'],
                                   coefficients[1:], policy['unmatched_label'])
    names = ('local', 'existence_bce', 'boundary_ce', 'termination_bce')
    losses = {'local': local_loss(outputs['segmentation_logits'][0], targets),
              **{name: query['components'][name] for name in names[1:]}}
    return {'total': weight*sum(c*losses[name] for c, name in zip(coefficients, names)),
            'losses': losses, 'query': query, **policy}
