"""A11 fixed-score TRAIN probe; reference-free, no training or final selector.

The weighted search is run once. Its sum/mean completed pools therefore see
identical histories; only completion retention differs. All three scores are
non-biological pointer/EOS receipts, not annotation gains or probabilities.
"""
import math
import torch
from torch.nn import functional as F
from .candidates import Budget, ORFIndex
from .query_decode import _end_events, _legal_ends
from .query_pointer import STAGES

ARMS = ('RAW_SUM', 'COMPONENT_SUM', 'COMPONENT_MEAN')


def completed_score(row, arm):
    if arm == 'RAW_SUM':
        return row['path_log_score']
    weighted = row['boundary_log_sum'] + 2 * row['termination_log_sum']
    if arm == 'COMPONENT_SUM':
        return weighted
    if arm == 'COMPONENT_MEAN' and row['chain']:
        return weighted / (2 * len(row['chain']))
    raise ValueError('Explicit complete-chain score arm required')


def _best(rows, limit, score):
    return sorted(rows, key=lambda r: (-score(r), r['chain'], r.get('entry', -1)))[:limit]


def _search(window, index, starts, acceptors, events, orf, limits, weighted):
    max_exons, beam_width, completed_limit = limits
    model = window.model
    arms = ARMS[1:] if weighted else ARMS[:1]
    pools = {arm: {} for arm in arms}
    counts = dict(entry_candidates=0, entry_pruned=0, legal_end_events=0,
                  completed_considered=0, continuation_pruned=0,
                  live_at_exon_cap=0, score_bound_stopped_live=0, depths=0)
    live = [dict(state=window.initial[index], chain=(), carry='', previous_end=-1,
                 path_log_score=0., boundary_log_sum=0., termination_log_sum=0.)]
    search_score = (lambda r: r['boundary_log_sum'] + 2*r['termination_log_sum']) if weighted else (
        lambda r: r['path_log_score'])
    for depth in range(max_exons):
        entries, entry_count = [], 0
        for row in live:
            lower = 0 if depth == 0 else row['previous_end'] + Budget().min_intron_bp
            coordinates = starts if depth == 0 else [a for a in acceptors if a >= lower]
            if not coordinates:
                continue
            stage = 'start' if depth == 0 else 'next_start'
            raw = F.log_softmax(window.pointer_logits(row['state'], stage, lower), dim=-1)
            scores = raw[torch.tensor(coordinates, device=raw.device)]
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal entry score')
            order = torch.argsort(scores, descending=True, stable=True)[:beam_width]
            entry_count += len(coordinates)
            for k in order.tolist():
                a, value = coordinates[k], float(scores[k])
                entries.append(dict(row, entry=a, state=window.advance(row['state'], stage, a),
                                    path_log_score=row['path_log_score']+value,
                                    boundary_log_sum=row['boundary_log_sum']+value))
        kept = _best(entries, beam_width, search_score)
        counts['entry_candidates'] += entry_count
        counts['entry_pruned'] += entry_count-len(kept)
        continuations = []
        for row in kept:
            a = row['entry']
            choices = _legal_ends(orf, events, a, row['carry'], first_exon=depth == 0)
            counts['legal_end_events'] += len(choices)
            if not choices:
                continue
            coordinates = torch.tensor([b for b, _, _ in choices], device=row['state'].device)
            tokens = (window.end_features[coordinates]+window.positions[coordinates]+
                      model.event.weight[STAGES['end']])
            states = model.recurrent(tokens, row['state'].expand(len(choices), -1))
            eos = model.terminate(states).squeeze(-1)
            terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices], device=eos.device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
            boundary = F.log_softmax(window.pointer_logits(row['state'], 'end', a), dim=-1)[coordinates]
            combined = boundary+branch
            if not bool(torch.isfinite(combined).all()):
                raise ValueError('Nonfinite legal end/EOS score')
            for k, (b, kind, carry) in enumerate(choices):
                new = dict(state=states[k], chain=row['chain']+((a, b),), carry=carry,
                           previous_end=b, path_log_score=row['path_log_score']+float(combined[k]),
                           boundary_log_sum=row['boundary_log_sum']+float(boundary[k]),
                           termination_log_sum=row['termination_log_sum']+float(branch[k]))
                if kind == 'stop':
                    counts['completed_considered'] += 1
                    for arm in arms:
                        old = pools[arm].get(new['chain'])
                        if old is None or completed_score(new, arm) > completed_score(old, arm):
                            pools[arm][new['chain']] = new
                else:
                    continuations.append(new)
        for arm in arms:
            done = _best(list(pools[arm].values()), completed_limit,
                         lambda r: completed_score(r, arm))
            pools[arm] = {r['chain']: r for r in done}
        live = _best(continuations, beam_width, search_score)
        counts['continuation_pruned'] += len(continuations)-len(live)
        counts['depths'] = depth+1
        if not live:
            break
        if depth == max_exons-1:
            counts['live_at_exon_cap'] = len(live)
        # No score bound. Completion pools never affect the retained live beam.
    result = {}
    for arm in arms:
        rows = _best(list(pools[arm].values()), completed_limit, lambda r: completed_score(r, arm))
        completed = []
        for rank, row in enumerate(rows):
            receipt = {k: row[k] for k in ('chain', 'path_log_score', 'boundary_log_sum', 'termination_log_sum')}
            receipt.update(query_index=index, candidate_rank=rank,
                           query_logit=float(window.query_logits[index]),
                           completed_score=completed_score(row, arm), terminated=True,
                           stop_reason='canonical_stop')
            if not all(math.isfinite(receipt[k]) for k in
                       ('path_log_score', 'boundary_log_sum', 'termination_log_sum', 'completed_score')):
                raise ValueError('Nonfinite completion receipt')
            completed.append(receipt)
        result[arm] = dict(query_index=index, query_logit=float(window.query_logits[index]),
                           completed=completed, counts=dict(counts))
    return result


@torch.no_grad()
def decode_query_component_pair(model, features, sequence, max_exons, beam_width,
                                completed_limit, valid_bases=None):
    if min(max_exons, beam_width, completed_limit) < 1:
        raise ValueError('Positive explicit exon/beam/completion limits required')
    window = model.for_window(features, valid_bases)
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Oriented DNA and features must align')
    if not bool(torch.isfinite(window.query_logits).all()):
        raise ValueError('Nonfinite query existence logit')
    starts = [a for a in range(len(sequence)-2) if sequence[a:a+3] == 'ATG']
    acceptors = [a+2 for a in range(len(sequence)-2) if sequence[a:a+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    result = {arm: [] for arm in ARMS}
    for q in range(model.num_queries):
        for weighted in (False, True):
            searched = _search(window, q, starts, acceptors, events, orf,
                               (max_exons, beam_width, completed_limit), weighted)
            for arm, rows in searched.items():
                result[arm].append(rows)
    return result
