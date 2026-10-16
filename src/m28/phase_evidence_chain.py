"""A26 frozen local-state odds under candidate CDS roles, no reference inputs."""
import math
import torch
from torch import nn
from .junction_chain import JunctionChainHead


def candidate_roles(chain, length, device=None):
    """Oriented half-open complete-CDS roles in the existing ascending phase order."""
    if not chain:
        raise ValueError('Null is not a complete CDS chain')
    first, last = int(chain[0][0]), int(chain[-1][1])
    roles = torch.zeros(length, dtype=torch.long, device=device)
    previous, coding = -1, 0
    for ordinal, (start, end) in enumerate(chain):
        if not 0 <= start < end <= length or (ordinal and start < previous):
            raise ValueError('Expected oriented nonoverlapping half-open exons')
        roles[start:end] = 1+(coding+torch.arange(end-start, device=device)) % 3
        if ordinal:
            roles[previous:start] = 4+(coding-1) % 3
            roles[previous-1] = 7+(coding-1) % 3
            roles[start] = 10+coding % 3
        coding += end-start
        previous = end
    if chain[0][1]-first < 3 or last-chain[-1][0] < 3:
        raise ValueError('Current raw grammar requires contiguous start and stop codons')
    # Match original fine-label precedence, including short terminal exons.
    roles[first:first+3] = 13
    roles[last-3:last] = 14
    return roles


@torch.no_grad()
def phase_evidence(local_logits, chains, valid_bases=None):
    """14 mean log odds against class0, then 14 log1p role counts.

    log P(state)-log P(0) equals the raw logit difference, avoiding an epsilon
    or extra normalization. Missing roles have mean0/count0. No labels enter.
    """
    length = len(local_logits) if valid_bases is None else int(valid_bases)
    if not 0 < length <= len(local_logits) or local_logits.shape[1] != 15:
        raise ValueError('Expected valid oriented 15-class local logits')
    local_logits = local_logits[:length]
    rows = []
    for chain in chains:
        roles = candidate_roles(chain, length, local_logits.device)
        means, counts = [], []
        for state in range(1, 15):
            keep = roles == state
            count = int(keep.sum())
            means.append((local_logits[keep, state]-local_logits[keep, 0]).mean()
                         if count else local_logits.new_tensor(0.))
            counts.append(math.log1p(count))
        rows.append(torch.cat([torch.stack(means), local_logits.new_tensor(counts)]))
    return torch.stack(rows) if rows else local_logits.new_empty((0, 28))


class PhaseEvidenceChainHead(JunctionChainHead):
    def __init__(self, hidden=128):
        super().__init__(hidden)
        self.path_evidence = nn.Linear(28, 1, bias=False)
        nn.init.zeros_(self.path_evidence.weight)

    def forward(self, features, chains, evidence, valid_bases=None):
        if list(evidence.shape) != [len(chains), 28]:
            raise ValueError('One ordered phase-evidence row per original chain required')
        scores = super().forward(features, chains, valid_bases)
        term = self.path_evidence(evidence).squeeze(-1)
        recurrent = scores['nonadditive']
        scores['nonadditive'] = recurrent+term
        scores['logits'] = scores['additive']+scores['nonadditive']-scores['null']
        scores['recurrent_nonadditive'] = recurrent
        scores['phase_evidence_term'] = term
        return scores
