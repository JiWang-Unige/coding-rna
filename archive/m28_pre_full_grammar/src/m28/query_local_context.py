"""A17 local query-context pooling; pointers and collaboration remain global."""
import math
import torch
from torch import nn
from .core import SharedContext
from .query_pointer import QueryPointerWindow
from .query_collaboration import CollaborativeGeneQueryModel, CollaborativeQueryPointer


def context_intervals(length, queries):
    # Valid-base intervals, not padded width or target-dependent positions.
    # For a short supported contig, adjacent queries may share one base.
    return tuple((i*length//queries, max(i*length//queries+1, (i+1)*length//queries))
                 for i in range(queries))


def local_context_attention(model, features, positions):
    keys = model.context_key(features)+positions
    logits = model.queries @ keys.T / math.sqrt(model.hidden)
    intervals = context_intervals(len(features), model.num_queries)
    coordinate = torch.arange(len(features), device=features.device)
    allowed = torch.stack([(coordinate >= a) & (coordinate < b) for a, b in intervals])
    return torch.softmax(logits.masked_fill(~allowed, float('-inf')), dim=-1)


class LocalContextQueryWindow(QueryPointerWindow):
    def __init__(self, model, features, valid_bases):
        # Reuse the unmodified geometric pointer keys and recurrent interface.
        # The original global initial state is replaced, never used downstream.
        super().__init__(model, features, valid_bases)
        attention = local_context_attention(model, self.features, self.positions[:-1])
        context = attention @ model.context_value(self.features)
        initial = model.initialize(torch.cat([model.queries, context], dim=-1))
        q = model.query_norm(initial)[None]
        mixed, _ = model.query_attention(q, q, q, need_weights=False)
        self.initial = initial+mixed[0]
        self.query_logits = model.exists(self.initial).squeeze(-1)


class LocalContextQueryPointer(CollaborativeQueryPointer):
    def for_window(self, features, valid_bases=None):
        return LocalContextQueryWindow(self, features.detach(), valid_bases)


class LocalContextGeneQueryModel(CollaborativeGeneQueryModel):
    def __init__(self, feature_dim, num_queries, hidden=128):
        nn.Module.__init__(self)
        # Exactly the A3 constructor order, with no new parameter/random draw.
        self.shared = SharedContext(feature_dim, hidden)
        self.segmentation = nn.Linear(hidden, 15)
        self.query = LocalContextQueryPointer(num_queries, hidden)

    def query_window(self, hidden, valid_bases=None):
        return LocalContextQueryWindow(self.query, hidden, valid_bases)
