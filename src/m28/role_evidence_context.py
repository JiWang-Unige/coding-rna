"""Candidate-role evidence enters each exon/junction token before recurrence."""
import math
import torch
from torch import nn
from .core import SharedContext
from .phase_evidence_chain import candidate_roles


@torch.no_grad()
def role_margins(local_logits, chain, length):
    roles = candidate_roles(chain, length, local_logits.device)
    positions = torch.arange(length, device=local_logits.device)
    return local_logits[:length][positions, roles]-local_logits[:length, 0]


class RoleConditionedChainHead(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        # Original branch order, with 3 exon and 5 junction evidence columns.
        self.exon = nn.Sequential(nn.Linear(hidden*3+6, hidden), nn.GELU())
        self.sequence = nn.GRU(hidden, hidden, batch_first=True)
        self.additive = nn.Linear(hidden, 1)
        self.nonadditive = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh(), nn.Linear(hidden, 1))
        self.null = nn.Linear(hidden, 1)
        self.junction = nn.Sequential(nn.Linear(hidden*5+7, hidden), nn.GELU())
        self.junction_additive = nn.Linear(hidden, 1)

    def forward(self, features, chains, local_logits, valid_bases=None):
        length = len(features) if valid_bases is None else int(valid_bases)
        if (not 0 < length <= len(features) or local_logits.ndim != 2 or
                local_logits.shape[1] != 15 or len(local_logits) < length or
                local_logits.dtype != torch.float32 or not bool(torch.isfinite(local_logits[:length]).all())):
            raise ValueError('Valid original oriented FP32 15-class local evidence required')
        features = features[:length]
        null = self.null(features.mean(dim=0)).squeeze(-1)
        exon_add, junction_add, nonadd = [], [], []
        for chain in chains:
            margin = role_margins(local_logits, chain, length)
            exons, junctions = [], []
            previous = coding = 0
            for ordinal, (start, end) in enumerate(chain):
                gap = start-previous if ordinal else 0
                exon_geo = features.new_tensor([math.log1p(end-start), math.log1p(gap), coding % 3/2.])
                exons.append(torch.cat([features[start:end].mean(dim=0), features[start], features[end-1],
                    exon_geo, torch.stack([margin[start:end].mean(), margin[start], margin[end-1]])]))
                if ordinal:
                    zero = features.new_zeros(features.shape[1])
                    intron_mean = features[previous:start].mean(dim=0) if gap else zero
                    first = features[previous] if gap else zero
                    last = features[start-1] if gap else zero
                    junction_geo = features.new_tensor([math.log1p(gap), coding % 3/2.])
                    empty = margin.new_tensor(0.)
                    evidence = torch.stack([margin[previous:start].mean() if gap else empty,
                        margin[previous-1], margin[previous] if gap else empty,
                        margin[start-1] if gap else empty, margin[start]])
                    junctions.append(torch.cat([intron_mean, features[previous-1], first, last,
                        features[start], junction_geo, evidence]))
                previous = end; coding += end-start
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


class RoleEvidenceContext(nn.Module):
    def __init__(self, feature_dim=2048, hidden=128):
        super().__init__()
        self.shared = SharedContext(feature_dim, hidden)
        self.head = RoleConditionedChainHead(hidden)

    def forward(self, features, dna, valid_bases, chains, local_logits):
        if dna.shape[0] != 1 or len(valid_bases) != 1:
            raise ValueError('One original oriented window and its candidate pool required')
        hidden, _ = self.shared(features, dna, valid_bases)
        return self.head(hidden[0], chains, local_logits, valid_bases[0])
