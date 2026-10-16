"""A21 complete-chain score with explicit interleaved intron/junction tokens.

This changes the head representation, not candidate generation or annotation
semantics. Coordinates and features are in the existing strand-oriented frame.
"""
import math
import torch
from torch import nn


class JunctionChainHead(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        # Keep these five branches in the original fresh-seed initialization order.
        self.exon = nn.Sequential(nn.Linear(hidden*3+3, hidden), nn.GELU())
        self.sequence = nn.GRU(hidden, hidden, batch_first=True)
        self.additive = nn.Linear(hidden, 1)
        self.nonadditive = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh(), nn.Linear(hidden, 1))
        self.null = nn.Linear(hidden, 1)
        self.junction = nn.Sequential(nn.Linear(hidden*5+2, hidden), nn.GELU())
        self.junction_additive = nn.Linear(hidden, 1)

    def forward(self, features, chains, valid_bases=None):
        length = len(features) if valid_bases is None else int(valid_bases)
        if not 0 < length <= len(features):
            raise ValueError('Invalid valid length for oriented chain features')
        features = features[:length]
        null = self.null(features.mean(dim=0)).squeeze(-1)
        exon_add, junction_add, nonadd = [], [], []
        for chain in chains:
            if not chain:
                raise ValueError('The empty chain is null, not an exon list')
            exons, junctions = [], []
            previous, coding = 0, 0
            for ordinal, (start, end) in enumerate(chain):
                if not 0 <= start < end <= length or (ordinal and start < previous):
                    raise ValueError('Expected ordered nonoverlapping half-open exons')
                gap = start-previous if ordinal else 0
                exon_geo = features.new_tensor([math.log1p(end-start), math.log1p(gap), coding % 3/2.])
                exons.append(torch.cat([features[start:end].mean(dim=0),
                                        features[start], features[end-1], exon_geo]))
                if ordinal:
                    # Empty gaps are legal in the existing candidate contract.
                    # No intronic bases exist there: use zero intron features.
                    intron_mean = features[previous:start].mean(dim=0) if gap else features.new_zeros(features.shape[1])
                    first = features[previous] if gap else features.new_zeros(features.shape[1])
                    last = features[start-1] if gap else features.new_zeros(features.shape[1])
                    junction_geo = features.new_tensor([math.log1p(gap), coding % 3/2.])
                    junctions.append(torch.cat([intron_mean, features[previous-1],
                                                first, last, features[start], junction_geo]))
                previous = end
                coding += end-start
            exon_tokens = self.exon(torch.stack(exons))
            junction_tokens = self.junction(torch.stack(junctions)) if junctions else None
            tokens = []
            for ordinal, token in enumerate(exon_tokens):
                if ordinal:
                    tokens.append(junction_tokens[ordinal-1])
                tokens.append(token)
            _, state = self.sequence(torch.stack(tokens)[None])
            exon_add.append(self.additive(exon_tokens).sum())
            junction_add.append(self.junction_additive(junction_tokens).sum()
                                if junctions else exon_tokens.sum()*0)
            nonadd.append(self.nonadditive(state[-1, 0]).squeeze(-1))
        exon = torch.stack(exon_add) if exon_add else features.new_empty((0,))
        junction = torch.stack(junction_add) if junction_add else features.new_empty((0,))
        non = torch.stack(nonadd) if nonadd else features.new_empty((0,))
        additive = exon+junction
        return dict(logits=additive+non-null, additive=additive, nonadditive=non,
                    null=null, additive_only_logits=additive-null,
                    exon_additive=exon, junction_additive=junction)
