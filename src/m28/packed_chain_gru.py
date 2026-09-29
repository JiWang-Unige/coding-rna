"""Isolated chain-GRU batching prototype; not imported by the R4 runners.

Only the GRU call is batched. Pooling, exon projection, readouts, parameters,
candidate order and valid-base handling follow core.ChainHead.forward.
"""
import math
import torch
from torch.nn.utils.rnn import pack_sequence


def packed_chain_forward(head, features, chains, valid_bases=None):
    length = len(features) if valid_bases is None else int(valid_bases)
    features = features[:length]
    null = head.null(features.mean(dim=0)).squeeze(-1)
    sequences, additive = [], []
    for chain in chains:
        if not chain:
            raise ValueError("Empty chain is represented by null, not by an exon list")
        values = []
        previous = coding = 0
        for j, (a, b) in enumerate(chain):
            if not 0 <= a < b <= length or (j and a < previous):
                raise ValueError("Exons must be nonoverlapping, ordered, oriented half-open intervals")
            geometry = features.new_tensor([
                math.log1p(b-a), math.log1p(a-previous) if j else 0., coding % 3 / 2.])
            values.append(torch.cat([
                features[a:b].mean(dim=0), features[a], features[b-1], geometry]))
            previous = b
            coding += b-a
        z = head.exon(torch.stack(values))
        sequences.append(z)
        additive.append(head.additive(z).sum())
    if sequences:
        # Each sequence has its own implicit zero initial state; packing excludes
        # padding. GRU returns h in original input order for an unsorted pack.
        _, h = head.sequence(pack_sequence(sequences, enforce_sorted=False))
        nonadditive = [
            head.nonadditive(h[-1, i]).squeeze(-1) for i in range(len(sequences))]
        add = torch.stack(additive)
        non = torch.stack(nonadditive)
    else:
        add = features.new_empty((0,))
        non = features.new_empty((0,))
    return {"additive": add, "nonadditive": non, "null": null,
            "logits": add+non-null, "additive_only_logits": add-null}
