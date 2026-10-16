"""A37 fresh context learned directly through ordered complete-chain evidence."""
from torch import nn
from .core import SharedContext
from .junction_chain import JunctionChainHead


class ContextChain(nn.Module):
    def __init__(self, feature_dim=2048, hidden=128):
        super().__init__()
        self.shared = SharedContext(feature_dim, hidden)
        self.head = JunctionChainHead(hidden)

    def forward(self, features, dna, valid_bases, chains):
        if dna.shape[0] != 1 or len(valid_bases) != 1:
            raise ValueError('One original oriented window and its candidate pool required')
        hidden, _ = self.shared(features, dna, valid_bases)
        return self.head(hidden[0], chains, valid_bases[0])
