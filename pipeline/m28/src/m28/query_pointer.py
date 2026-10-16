"""A-backup gene queries and boundary pointers, synthetic interface only.

This raw geometric interface is not a FASTA annotation caller. Canonical DNA/
ORF constraints are provided separately in query_decode; one-to-one assignment
primitives live in query_training. Real-data matching/label policies, loss
weights and the first-fit protocol remain unset. Known-query labels are supplied
explicitly; missing queries are never automatically labelled negative.
Window state is fresh per optimizer step.
"""
import math
import torch
from torch import nn

STAGES = {'start': 0, 'end': 1, 'next_start': 2}


class GeneQueryPointer(nn.Module):
    def __init__(self, num_queries, hidden=128):
        super().__init__()
        self.num_queries, self.hidden = num_queries, hidden
        self.queries = nn.Parameter(torch.empty(num_queries, hidden))
        nn.init.normal_(self.queries, std=0.02)
        self.context_key = nn.Linear(hidden, hidden)
        self.context_value = nn.Linear(hidden, hidden)
        self.coordinate = nn.Linear(1, hidden, bias=False)
        self.initialize = nn.Sequential(nn.Linear(2*hidden, hidden), nn.Tanh())
        self.exists = nn.Linear(hidden, 1)
        self.start_key = nn.Linear(hidden, hidden, bias=False)
        self.end_key = nn.Linear(hidden, hidden, bias=False)
        self.event = nn.Embedding(3, hidden)
        self.pointer_query = nn.Linear(hidden, hidden, bias=False)
        self.recurrent = nn.GRUCell(hidden, hidden)
        self.terminate = nn.Linear(hidden, 1)

    def for_window(self, features, valid_bases=None):
        # The raw interface continues to freeze externally supplied features.
        return QueryPointerWindow(self, features.detach(), valid_bases)

    def teacher_components(self, features, query_labels, positive_chains, valid_bases=None):
        return self.for_window(features, valid_bases).teacher_components(query_labels, positive_chains)

    @torch.no_grad()
    def greedy_queries(self, features, max_exons, valid_bases=None):
        window = self.for_window(features, valid_bases)
        return [window.greedy_query(i, max_exons) for i in range(self.num_queries)]


class QueryPointerWindow:
    def __init__(self, model, features, valid_bases):
        self.model = model
        self.length = len(features) if valid_bases is None else int(valid_bases)
        if not 0 < self.length <= len(features):
            raise ValueError('Expected a nonempty valid window inside its feature tensor')
        # A's new shared-context caller retains its own trainable graph here.
        self.features = features[:self.length]
        position = torch.arange(self.length+1, device=features.device, dtype=features.dtype)[:, None]/self.length
        self.positions = model.coordinate(position)
        keys = model.context_key(self.features) + self.positions[:-1]
        attention = torch.softmax(model.queries @ keys.T / math.sqrt(model.hidden), dim=-1)
        context = attention @ model.context_value(self.features)
        self.initial = model.initialize(torch.cat([model.queries, context], dim=-1))
        self.query_logits = model.exists(self.initial).squeeze(-1)
        # All coordinates are half-open boundaries 0..L. Entry L is invalid;
        # exon end L is valid and uses the preceding base's feature.
        self.entry_features = torch.cat([self.features, self.features[-1:]])
        self.end_features = torch.cat([self.features[:1], self.features])
        self.entry_keys = model.start_key(self.entry_features) + self.positions
        self.end_keys = model.end_key(self.end_features) + self.positions

    def pointer_logits(self, state, stage, lower):
        event = self.model.event.weight[STAGES[stage]]
        query = self.model.pointer_query(state + event)
        keys = self.end_keys if stage == 'end' else self.entry_keys
        # Keep the original scalar path for teacher/greedy decoding.
        logits = ((keys @ query) if state.ndim == 1 else (query @ keys.T))
        logits = logits / math.sqrt(self.model.hidden)
        coordinate = torch.arange(self.length+1, device=logits.device)
        allowed = ((coordinate > lower) & (coordinate <= self.length) if stage == 'end'
                   else (coordinate >= lower) & (coordinate < self.length))
        return logits.masked_fill(~allowed, float('-inf'))

    def advance(self, state, stage, coordinate):
        local = self.end_features[coordinate] if stage == 'end' else self.entry_features[coordinate]
        token = local + self.positions[coordinate] + self.model.event.weight[STAGES[stage]]
        if state.ndim == 2:
            token = token.expand_as(state)
        return self.model.recurrent(token, state)

    def validate_chain(self, chain):
        if not chain:
            raise ValueError('A positive query requires a nonempty exon chain')
        previous = -1
        for a, b in chain:
            if not 0 <= a < b <= self.length or a <= previous:
                raise ValueError('Expected ordered oriented half-open exons with positive gaps')
            previous = b

    def teacher_components(self, query_labels, positive_chains):
        labels = torch.as_tensor(query_labels, device=self.features.device, dtype=self.features.dtype)
        if labels.shape != (self.model.num_queries,) or not ((labels == -1)|(labels == 0)|(labels == 1)).all():
            raise ValueError('One explicit known-positive, known-negative or unknown label per query is required')
        expected = {i for i, value in enumerate(labels.tolist()) if value == 1}
        if set(positive_chains) != expected:
            raise ValueError('Positive query assignment must match the supplied query labels')
        zero = self.initial.sum()*0 + self.query_logits.sum()*0
        known = labels >= 0
        existence = (nn.functional.binary_cross_entropy_with_logits(self.query_logits[known], labels[known])
                     if known.any() else zero)
        boundary_losses, termination_losses = [], []
        for i in sorted(positive_chains):
            chain = positive_chains[i]
            self.validate_chain(chain)
            state, previous_end = self.initial[i], -1
            boundaries, endings = [], []
            for j, (a, b) in enumerate(chain):
                stage = 'start' if j == 0 else 'next_start'
                lower = 0 if j == 0 else previous_end+1
                boundaries.append(nn.functional.cross_entropy(
                    self.pointer_logits(state, stage, lower)[None],
                    torch.tensor([a], device=labels.device)))
                state = self.advance(state, stage, a)
                boundaries.append(nn.functional.cross_entropy(
                    self.pointer_logits(state, 'end', a)[None],
                    torch.tensor([b], device=labels.device)))
                state = self.advance(state, 'end', b)
                endings.append(nn.functional.binary_cross_entropy_with_logits(
                    self.model.terminate(state).squeeze(-1), labels.new_tensor(float(j == len(chain)-1))))
                previous_end = b
            boundary_losses.append(torch.stack(boundaries).mean())
            termination_losses.append(torch.stack(endings).mean())
        # Separate components: no first-fit coefficients are silently chosen here.
        return {'existence_bce': existence,
                'boundary_ce': torch.stack(boundary_losses).mean() if boundary_losses else zero,
                'termination_bce': torch.stack(termination_losses).mean() if termination_losses else zero,
                'known_queries': int(known.sum()), 'positive_queries': len(positive_chains)}

    def greedy_query(self, index, max_exons):
        state, previous_end, chain = self.initial[index], -1, []
        def result(terminated, reason):
            return {'query_index': index, 'query_logit': float(self.query_logits[index]),
                    'chain': tuple(chain), 'terminated': terminated, 'stop_reason': reason}
        for j in range(max_exons):
            lower = 0 if j == 0 else previous_end+1
            if lower >= self.length:
                return result(False, 'no_next_start')
            stage = 'start' if j == 0 else 'next_start'
            a = int(self.pointer_logits(state, stage, lower).argmax())
            state = self.advance(state, stage, a)
            b = int(self.pointer_logits(state, 'end', a).argmax())
            state = self.advance(state, 'end', b)
            chain.append((a, b))
            if float(self.model.terminate(state)) > 0:
                return result(True, 'predicted_termination')
            previous_end = b
        return result(False, 'exon_limit')
