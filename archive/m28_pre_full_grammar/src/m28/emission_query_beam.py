"""A33 global beam retention under the fixed A31 emission objective.

The original local pointer top-k and DNA grammar remain unchanged. This is
bounded search, not an exact decoder, new training or a calibrated likelihood.
"""
import math
import torch
from torch.nn import functional as F
from .candidates import Budget, ORFIndex
from .query_decode import _end_events, _legal_ends
from .query_pointer import STAGES


class EmissionScores:
    """FP64 interval sums with the original start/stop/donor/acceptor precedence."""
    def __init__(self, logits, length):
        if logits.ndim != 2 or logits.shape[1] != 15 or not 0 < length <= len(logits):
            raise ValueError('Valid oriented 15-class emissions required')
        gain = (logits[:length].double()-logits[:length, :1].double()).cpu()
        if not bool(torch.isfinite(gain).all()):
            raise ValueError('Nonfinite valid emission')
        self.length = length
        self.gain = gain.tolist()
        self.prefix = torch.cat((torch.zeros(1, 15, dtype=torch.float64), gain.cumsum(0))).T.tolist()
        positions = torch.arange(length)
        coding = torch.stack([gain[positions, 1+(positions+offset) % 3] for offset in range(3)], dim=1)
        self.coding = torch.cat((torch.zeros(1, 3, dtype=torch.float64), coding.cumsum(0))).T.tolist()
        positive = gain[:, 1:].max(1).values.clamp_min(0.)
        self.suffix = torch.cat((positive.flip(0).cumsum(0).flip(0), torch.zeros(1, dtype=torch.float64))).tolist()

    def interval(self, a, b, role):
        return self.prefix[role][b]-self.prefix[role][a]

    def coding_interval(self, a, b, coding):
        offset = (coding-a) % 3
        return self.coding[offset][b]-self.coding[offset][a]

    def entry_gain(self, row, a):
        if not row['chain']:
            return self.interval(a, a+3, 13)
        phase = row['coding_bases'] % 3
        return (row['emission_gain']+self.interval(row['previous_end'], a, 4+(phase-1) % 3)+
                self.gain[a][10+phase])

    def extend(self, row, a, b, terminal):
        size, coding = b-a, row['coding_bases']
        first = not row['chain']
        value = row['emission_gain']+self.coding_interval(a, b, coding)
        if first:
            value += self.interval(a, a+3, 13)-self.coding_interval(a, a+3, 0)
        else:
            value += self.interval(row['previous_end'], a, 4+(coding-1) % 3)
            # A contiguous three-base terminal STOP overrides the acceptor.
            if not (terminal and size == 3):
                value += self.gain[a][10+coding % 3]-self.gain[a][1+coding % 3]
        if terminal:
            value += self.interval(b-3, b, 14)-self.coding_interval(b-3, b, coding+size-3)
            # Only an artificial three-base single exon overlaps START/STOP.
            if first and size == 3:
                value -= self.interval(a, b, 13)-self.coding_interval(a, b, 0)
        elif not (first and size == 3):
            phase = (coding+size-1) % 3
            prior = 10+coding % 3 if not first and size == 1 else 1+phase
            value += self.gain[b-1][7+phase]-self.gain[b-1][prior]
        return value

    def upper_bound(self, row):
        # Every unseen position contributes at most its largest positive role
        # gain. Only retained histories are bounded; earlier pruning is not undone.
        return row['emission_gain']+self.suffix[row['previous_end']]


def _best(rows, limit):
    return sorted(rows, key=lambda r: (-r['search_gain'], -r['path_log_score'],
                                      r['chain'], r.get('entry', -1)))[:limit]


def _emission_query(window, index, starts, acceptors, events, orf, scorer,
                    max_exons, beam_width, completed_limit):
    model = window.model
    counts = dict(entry_candidates=0, entry_pruned=0, legal_end_events=0,
        completed_considered=0, continuation_pruned=0, live_at_exon_cap=0,
        emission_bound_stopped_live=0, depths_visited=0)
    live = [dict(state=window.initial[index], chain=(), carry='', previous_end=-1,
                 coding_bases=0, emission_gain=0., search_gain=0., path_log_score=0.)]
    completed = {}
    for depth in range(max_exons):
        counts['depths_visited'] += 1
        entries = []; considered = 0
        for row in live:
            lower = 0 if depth == 0 else row['previous_end']+Budget().min_intron_bp
            coordinates = starts if depth == 0 else [a for a in acceptors if a >= lower]
            if not coordinates:
                continue
            stage = 'start' if depth == 0 else 'next_start'
            logprob = F.log_softmax(window.pointer_logits(row['state'], stage, lower), dim=-1)
            values = logprob[torch.tensor(coordinates, device=logprob.device)]
            if not bool(torch.isfinite(values).all()):
                raise ValueError('Nonfinite legal entry score')
            # Local geometric pointer top-k is deliberately unchanged.
            order = torch.argsort(values, descending=True, stable=True)[:beam_width]
            counts['entry_candidates'] += len(coordinates); considered += len(coordinates)
            for k in order.tolist():
                a = coordinates[k]
                entries.append({**row, 'entry': a, 'state': window.advance(row['state'], stage, a),
                    'path_log_score': row['path_log_score']+float(values[k]),
                    'search_gain': scorer.entry_gain(row, a)})
        kept = _best(entries, beam_width)
        counts['entry_pruned'] += considered-len(kept)
        continuations = []
        for row in kept:
            a = row['entry']
            choices = _legal_ends(orf, events, a, row['carry'], first_exon=depth == 0)
            counts['legal_end_events'] += len(choices)
            if not choices:
                continue
            coordinates = torch.tensor([b for b, _, _ in choices], device=row['state'].device)
            tokens = window.end_features[coordinates]+window.positions[coordinates]+model.event.weight[STAGES['end']]
            states = model.recurrent(tokens, row['state'].expand(len(choices), -1))
            eos = model.terminate(states).squeeze(-1)
            terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices], device=row['state'].device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
            values = F.log_softmax(window.pointer_logits(row['state'], 'end', a), dim=-1)[coordinates]+branch
            if not bool(torch.isfinite(values).all()):
                raise ValueError('Nonfinite legal end/EOS score')
            for k, (b, kind, carry) in enumerate(choices):
                gain = scorer.extend(row, a, b, kind == 'stop')
                if not math.isfinite(gain):
                    raise ValueError('Nonfinite prefix emission gain')
                new = dict(state=states[k], chain=row['chain']+((a, b),), carry=carry,
                    previous_end=b, coding_bases=row['coding_bases']+b-a, emission_gain=gain,
                    search_gain=gain, path_log_score=row['path_log_score']+float(values[k]))
                if kind == 'stop':
                    counts['completed_considered'] += 1
                    old = completed.get(new['chain'])
                    if old is None or new['path_log_score'] > old['path_log_score']:
                        completed[new['chain']] = new
                else:
                    continuations.append(new)
        done = _best(list(completed.values()), completed_limit)
        completed = {r['chain']: r for r in done}
        live = _best(continuations, beam_width)
        counts['continuation_pruned'] += len(continuations)-len(live)
        if not live:
            break
        if depth == max_exons-1:
            counts['live_at_exon_cap'] = len(live); break
        # Do not reuse the old nonpositive log-increment stopping rule: future
        # role gains can be positive. Strict inequality retains emission ties.
        if len(done) == completed_limit and max(scorer.upper_bound(r) for r in live) < done[-1]['emission_gain']:
            counts['emission_bound_stopped_live'] = len(live); break
    results = [dict(query_index=index, candidate_rank=rank, query_logit=float(window.query_logits[index]),
        chain=row['chain'], path_log_score=row['path_log_score'], emission_gain=row['emission_gain'],
        terminated=True, stop_reason='canonical_stop')
        for rank, row in enumerate(_best(list(completed.values()), completed_limit))]
    return dict(query_index=index, query_logit=float(window.query_logits[index]), completed=results, counts=counts)


@torch.no_grad()
def decode_emission_query_beam(model, features, logits, sequence, max_exons, beam_width,
                               completed_limit, valid_bases=None):
    if min(max_exons, beam_width, completed_limit) < 1:
        raise ValueError('Positive explicit search limits required')
    window = model.for_window(features, valid_bases)
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Oriented DNA and valid features differ')
    if not bool(torch.isfinite(window.query_logits).all()):
        raise ValueError('Nonfinite query existence logit')
    scorer = EmissionScores(logits, len(sequence))
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    return [_emission_query(window, i, starts, acceptors, events, orf, scorer,
            max_exons, beam_width, completed_limit) for i in range(model.num_queries)]


def augment_pool(chains, sources, queries):
    origins = {tuple(map(tuple, c)): set(s) for c, s in zip(chains, sources)}
    if len(chains) != len(sources) or len(origins) != len(chains):
        raise ValueError('Unique original geometry/source pairs required')
    for query in queries:
        for row in query['completed']:
            origins.setdefault(tuple(map(tuple, row['chain'])), set()).add('emission_beam4')
    ordered = sorted(origins)
    return [list(chain) for chain in ordered], [sorted(origins[chain]) for chain in ordered]
