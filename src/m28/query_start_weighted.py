"""A25 initiation prior in A17 local attention; chain pointers stay global."""
import math
import torch
from torch import nn
from .core import SharedContext
from .query_pointer import QueryPointerWindow
from .query_local_context import (context_intervals, LocalContextGeneQueryModel,
                                  LocalContextQueryPointer)
from .query_known_mean import joint_known_query_losses


def start_weighted_attention(model, features, positions, initiation_logits):
    keys = model.context_key(features)+positions
    logits = model.queries @ keys.T / math.sqrt(model.hidden)
    coordinate = torch.arange(len(features), device=features.device)
    allowed = torch.stack([(coordinate >= a) & (coordinate < b)
                           for a, b in context_intervals(len(features), model.num_queries)])
    return torch.softmax((logits+initiation_logits[None]).masked_fill(~allowed, float('-inf')), dim=-1)


class StartWeightedQueryWindow(QueryPointerWindow):
    def __init__(self, model, features, valid_bases):
        super().__init__(model, features, valid_bases)
        self.initiation_logits = model.initiation(self.features).squeeze(-1)
        attention = start_weighted_attention(model, self.features, self.positions[:-1], self.initiation_logits)
        context = attention @ model.context_value(self.features)
        initial = model.initialize(torch.cat([model.queries, context], dim=-1))
        q = model.query_norm(initial)[None]
        mixed, _ = model.query_attention(q, q, q, need_weights=False)
        self.initial = initial+mixed[0]
        self.query_logits = model.exists(self.initial).squeeze(-1)


class StartWeightedQueryPointer(LocalContextQueryPointer):
    def __init__(self, num_queries, hidden=128):
        super().__init__(num_queries, hidden)
        self.initiation = nn.Linear(hidden, 1, bias=False)
        nn.init.zeros_(self.initiation.weight)

    def for_window(self, features, valid_bases=None):
        return StartWeightedQueryWindow(self, features.detach(), valid_bases)


class StartWeightedGeneQueryModel(LocalContextGeneQueryModel):
    def __init__(self, feature_dim, num_queries, hidden=128):
        nn.Module.__init__(self)
        self.shared = SharedContext(feature_dim, hidden)
        self.segmentation = nn.Linear(hidden, 15)
        self.query = StartWeightedQueryPointer(num_queries, hidden)

    def query_window(self, hidden, valid_bases=None):
        return StartWeightedQueryWindow(self.query, hidden, valid_bases)


def initiation_components(logits, targets):
    """Conditional ranking of eligible first CDS bases, not gene absence BCE.

    Other fine-class-13 bases are not negatives. Fine unknown/padding positions
    are excluded unless an explicit eligible positive chain supplies that start.
    These masks affect this loss only, never the reference-free attention.
    """
    starts = tuple(sorted({int(chain[0][0]) for chain in targets['positive_chains']}))
    fine = torch.as_tensor(targets['seg'][:len(logits)], device=logits.device)
    support = (fine >= 0) & (fine != 13)
    if not starts:
        return dict(loss=logits.sum()*0, starts=starts, support=support)
    indices = torch.tensor(starts, device=logits.device)
    support[indices] = True
    loss = -logits.masked_fill(~support, float('-inf')).log_softmax(0)[indices].mean()
    return dict(loss=loss, starts=starts, support=support)


def joint_start_weighted_losses(model, outputs, targets, valid_bases, window_weight, coefficients):
    result = joint_known_query_losses(model, outputs, targets, valid_bases, window_weight, coefficients)
    logits = model.query.initiation(outputs['features'][0, :valid_bases]).squeeze(-1)
    start = initiation_components(logits, targets)
    result['losses']['initiation_ce'] = start['loss']
    result['total'] = result['total']+float(window_weight)*start['loss']
    result['initiation_starts'] = start['starts']
    result['initiation_support_bases'] = int(start['support'].sum())
    return result
