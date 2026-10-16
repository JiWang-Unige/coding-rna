"""A23 relative boundary prior, without changing A3 starts or DNA decoding."""
import math
import torch
from torch import nn
from .core import SharedContext
from .query_collaboration import (CollaborativeGeneQueryModel, CollaborativeQueryPointer,
                                  CollaborativeQueryWindow)


def relative_geometry(window, lower):
    distance = (torch.arange(window.length+1, device=window.features.device,
                             dtype=window.features.dtype)-lower).clamp_min(0)
    return torch.stack((distance/window.length,
                        torch.log1p(distance)/math.log1p(window.length)), dim=-1)


class RelativeBoundaryQueryWindow(CollaborativeQueryWindow):
    def pointer_logits(self, state, stage, lower):
        logits = super().pointer_logits(state, stage, lower)
        if stage == 'start':
            return logits
        column = 0 if stage == 'end' else 1
        return logits+self.model.relative_boundary(relative_geometry(self, lower))[:, column]


class RelativeBoundaryQueryPointer(CollaborativeQueryPointer):
    def __init__(self, num_queries, hidden=128):
        super().__init__(num_queries, hidden)
        self.relative_boundary = nn.Sequential(nn.Linear(2, 32), nn.GELU(), nn.Linear(32, 2, bias=False))
        # All common A3 parameters precede new random draws. A constant output
        # bias would cancel in the pointer softmax and is deliberately absent.
        nn.init.zeros_(self.relative_boundary[-1].weight)

    def for_window(self, features, valid_bases=None):
        return RelativeBoundaryQueryWindow(self, features.detach(), valid_bases)


class RelativeBoundaryGeneQueryModel(CollaborativeGeneQueryModel):
    def __init__(self, feature_dim, num_queries, hidden=128):
        nn.Module.__init__(self)
        self.shared = SharedContext(feature_dim, hidden)
        self.segmentation = nn.Linear(hidden, 15)
        self.query = RelativeBoundaryQueryPointer(num_queries, hidden)

    def query_window(self, hidden, valid_bases=None):
        return RelativeBoundaryQueryWindow(self.query, hidden, valid_bases)
