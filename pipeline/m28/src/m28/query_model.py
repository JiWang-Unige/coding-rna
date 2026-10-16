"""A trainable shared context and class-balanced query supervision.

External cached GLM/DNA tensors are fixed inputs; query losses reach the new
shared module. This is an untrained single-window implementation of the pending
first-fit proposal, not a genome caller or authorization for real fitting.
"""
import math
import torch
from torch import nn
from torch.nn import functional as F
from scipy.optimize import linear_sum_assignment
from .core import SharedContext
from .labels import local_loss
from .query_pointer import GeneQueryPointer, QueryPointerWindow
from .query_training import teacher_chain_components


class GeneQueryModel(nn.Module):
    def __init__(self, feature_dim, num_queries, hidden=128):
        super().__init__()
        # Same common initialization order as C0, before extra query parameters.
        self.shared = SharedContext(feature_dim, hidden)
        self.segmentation = nn.Linear(hidden, 15)
        self.query = GeneQueryPointer(num_queries, hidden)

    def forward(self, token_features, dna, valid_bases):
        hidden, mask = self.shared(token_features.detach(), dna.detach(), valid_bases)
        return {'features': hidden, 'mask': mask,
                'segmentation_logits': self.segmentation(hidden)}

    def query_window(self, hidden, valid_bases=None):
        # Hidden belongs to this trainable shared context, not an external cache.
        return QueryPointerWindow(self.query, hidden, valid_bases)


def _balanced_existence(window, labels):
    labels = torch.as_tensor(labels, device=window.query_logits.device)
    parts = []
    for value in (0, 1):
        keep = labels == value
        if keep.any():
            parts.append(F.binary_cross_entropy_with_logits(
                window.query_logits[keep],
                window.query_logits.new_full((int(keep.sum()),), float(value))))
    return (torch.stack(parts).mean() if parts
            else window.initial.sum()*0 + window.query_logits.sum()*0)


def balanced_query_components(window, positive_chains, coefficients, unmatched_label):
    """Discrete matching minimizes the exact class-balanced window objective.

    Coefficients explicitly order existence, mean boundary and mean EOS.
    Unmatched query labels are caller-supplied 0 or -1; no partial/alternative
    targets are inferred here. All supplied positive chains must be assigned.
    """
    coefficients = tuple(float(c) for c in coefficients)
    if len(coefficients) != 3 or any(not math.isfinite(c) or c < 0 for c in coefficients):
        raise ValueError('Three finite nonnegative loss coefficients are required')
    if unmatched_label not in (-1, 0):
        raise ValueError('Unmatched queries require explicit unknown or negative labels')
    chains = tuple(tuple(tuple(exon) for exon in chain) for chain in positive_chains)
    n, q = len(chains), window.model.num_queries
    if n > q:
        raise ValueError('Query capacity cannot drop supplied positive chains')
    with torch.no_grad():
        costs = window.query_logits.new_empty((0, q))
        gene_to_query = ()
        if n:
            both_classes = unmatched_label == 0 and n < q
            positive_weight = (0.5 if both_classes else 1.)/n
            delta = positive_weight*F.softplus(-window.query_logits)
            if both_classes:
                delta = delta - 0.5*F.softplus(window.query_logits)/(q-n)
            rows = []
            for chain in chains:
                pair = teacher_chain_components(window, chain)
                rows.append(coefficients[0]*delta +
                            coefficients[1]*pair['boundary_ce']/n +
                            coefficients[2]*pair['termination_bce']/n)
            costs = torch.stack(rows)
            if not torch.isfinite(costs).all():
                raise ValueError('Nonfinite balanced query assignment cost')
            gene_indices, query_indices = linear_sum_assignment(costs.cpu().numpy())
            assignment = [None]*n
            for gene, query in zip(gene_indices, query_indices):
                assignment[int(gene)] = int(query)
            gene_to_query = tuple(assignment)
    labels = [unmatched_label]*q
    assigned = {}
    for gene, query in enumerate(gene_to_query):
        labels[query] = 1
        assigned[query] = chains[gene]
    # The old teacher still supplies geometry/event means. Only the explicit
    # class-balanced existence readout changes, matching the costs above.
    components = window.teacher_components(labels, assigned)
    components['existence_bce'] = _balanced_existence(window, labels)
    names = ('existence_bce', 'boundary_ce', 'termination_bce')
    return {'total': sum(c*components[k] for c, k in zip(coefficients, names)),
            'components': components, 'query_labels': tuple(labels),
            'positive_chains': assigned, 'gene_to_query': gene_to_query,
            'assignment_costs': costs}


def unmatched_query_policy(targets, valid_bases):
    """Pending A policy: protect an entire uncertain window, retaining positives."""
    n = int(valid_bases)
    if not 0 < n <= len(targets['callable_bases']):
        raise ValueError('Expected valid bases inside the supplied callable mask')
    positives = {tuple(tuple(exon) for exon in chain) for chain in targets['positive_chains']}
    protected = {tuple(tuple(exon) for exon in chain) for chain in targets['protected_chains']}
    reasons = []
    if targets['unknown_chain_spans']:
        reasons.append('unknown_chain_span')
    if protected - positives:
        reasons.append('complete_protected_nonpositive_chain')
    if not targets['callable_bases'][:n].all():
        reasons.append('noncallable_valid_base')
    return {'unmatched_label': -1 if reasons else 0, 'protection_reasons': tuple(reasons)}


def joint_query_losses(model, outputs, targets, valid_bases, window_weight, coefficients):
    """Single-window local plus class-balanced query loss, no optimizer or search."""
    coefficients = tuple(float(c) for c in coefficients)
    weight = float(window_weight)
    if len(coefficients) != 4 or any(not math.isfinite(c) or c < 0 for c in coefficients):
        raise ValueError('Four finite nonnegative joint coefficients are required')
    if not math.isfinite(weight) or weight < 0:
        raise ValueError('A finite nonnegative window weight is required')
    if outputs['features'].shape[0] != 1:
        raise ValueError('Expected a single-window query loss')
    window = model.query_window(outputs['features'][0], valid_bases)
    policy = unmatched_query_policy(targets, valid_bases)
    query = balanced_query_components(window, targets['positive_chains'],
                                      coefficients[1:], policy['unmatched_label'])
    names = ('local', 'existence_bce', 'boundary_ce', 'termination_bce')
    losses = {'local': local_loss(outputs['segmentation_logits'][0], targets),
              **{k: query['components'][k] for k in names[1:]}}
    return {'total': weight*sum(c*losses[k] for c, k in zip(coefficients, names)),
            'losses': losses, 'query': query, **policy}
