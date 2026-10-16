"""A2 query collaboration; old A constructors and checkpoint behavior are unchanged."""
import torch
from torch import nn
from .core import SharedContext
from .query_model import GeneQueryModel
from .query_pointer import GeneQueryPointer, QueryPointerWindow


class CollaborativeQueryWindow(QueryPointerWindow):
    def __init__(self, model, features, valid_bases):
        super().__init__(model, features, valid_bases)
        q = model.query_norm(self.initial)[None]
        mixed, _ = model.query_attention(q, q, q, need_weights=False)
        self.initial = self.initial + mixed[0]
        # Existence and all teacher/free pointer histories use this same state.
        self.query_logits = model.exists(self.initial).squeeze(-1)


class CollaborativeQueryPointer(GeneQueryPointer):
    def __init__(self, num_queries, hidden=128):
        # All old parameters are initialized before any new random draw.
        super().__init__(num_queries, hidden)
        self.query_norm = nn.LayerNorm(hidden)
        self.query_attention = nn.MultiheadAttention(hidden, 4, dropout=0., batch_first=True)

    def for_window(self, features, valid_bases=None):
        return CollaborativeQueryWindow(self, features.detach(), valid_bases)


class CollaborativeGeneQueryModel(GeneQueryModel):
    def __init__(self, feature_dim, num_queries, hidden=128):
        nn.Module.__init__(self)
        # Exact original A construction order; only the pointer adds new layers.
        self.shared = SharedContext(feature_dim, hidden)
        self.segmentation = nn.Linear(hidden, 15)
        self.query = CollaborativeQueryPointer(num_queries, hidden)

    def query_window(self, hidden, valid_bases=None):
        # Retain the trainable shared graph; the raw external-feature path detaches.
        return CollaborativeQueryWindow(self.query, hidden, valid_bases)
