"""Bounded reference-free branching for recurrent gene-query pointers.

Independent experimental decoder: does not modify frozen greedy callers.
Cumulative geometric-pointer/EOS log scores are not biological probabilities.
No training, length normalization, existence threshold or genome selector here.
"""
import math
import torch
from torch.nn import functional as F
from .candidates import Budget, ORFIndex
from .query_decode import _end_events, _legal_ends
from .query_pointer import STAGES


def _best(rows, limit):
    return sorted(rows, key=lambda r: (-r['path_log_score'], r['chain'],
                                      r.get('entry', -1)))[:limit]


def _beam_query(window, index, starts, acceptors, events, orf,
                max_exons, beam_width, completed_limit):
    model = window.model
    counts = dict(entry_candidates=0, entry_pruned=0, legal_end_events=0,
                  completed_considered=0, continuation_pruned=0,
                  live_at_exon_cap=0, score_bound_stopped_live=0)
    live = [{'state': window.initial[index], 'chain': (), 'carry': '',
             'previous_end': -1, 'path_log_score': 0.}]
    completed = {}
    for depth in range(max_exons):
        entries = []
        depth_candidates = 0
        for row in live:
            lower = 0 if depth == 0 else row['previous_end'] + Budget().min_intron_bp
            coordinates = starts if depth == 0 else [a for a in acceptors if a >= lower]
            if not coordinates:
                continue
            stage = 'start' if depth == 0 else 'next_start'
            # Retain the original geometric normalizer; do not renormalize
            # over DNA motifs or silently equate this to the teacher mean loss.
            logprob = F.log_softmax(window.pointer_logits(row['state'], stage, lower), dim=-1)
            coordinate_tensor = torch.tensor(coordinates, device=logprob.device)
            scores = logprob[coordinate_tensor]
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal entry score')
            order = torch.argsort(scores, descending=True, stable=True)[:beam_width]
            counts['entry_candidates'] += len(coordinates)
            depth_candidates += len(coordinates)
            for k in order.tolist():
                a = coordinates[k]
                entries.append({**row, 'entry': a,
                                'state': window.advance(row['state'], stage, a),
                                'path_log_score': row['path_log_score'] + float(scores[k])})
        kept_entries = _best(entries, beam_width)
        # Include both per-history top-k and global entry pruning in the count.
        counts['entry_pruned'] += depth_candidates - len(kept_entries)
        continuations = []
        for row in kept_entries:
            a = row['entry']
            choices = _legal_ends(orf, events, a, row['carry'], first_exon=depth == 0)
            counts['legal_end_events'] += len(choices)
            if not choices:
                continue
            coordinates = torch.tensor([b for b, _, _ in choices], device=row['state'].device)
            tokens = (window.end_features[coordinates] + window.positions[coordinates] +
                      model.event.weight[STAGES['end']])
            states = model.recurrent(tokens, row['state'].expand(len(choices), -1))
            eos = model.terminate(states).squeeze(-1)
            terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices],
                                    device=row['state'].device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
            boundary = F.log_softmax(window.pointer_logits(row['state'], 'end', a), dim=-1)
            scores = boundary[coordinates] + branch
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal end/EOS score')
            for k, (b, kind, carry) in enumerate(choices):
                new = {'state': states[k], 'chain': row['chain'] + ((a, b),),
                       'carry': carry, 'previous_end': b,
                       'path_log_score': row['path_log_score'] + float(scores[k])}
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
        counts['continuation_pruned'] += len(continuations) - len(live)
        if not live:
            break
        if depth == max_exons - 1:
            counts['live_at_exon_cap'] = len(live)
            break
        # Every future log increment is <=0. This bound is valid for the
        # remaining retained beam, not a global optimum over pruned histories.
        if len(done) == completed_limit and live[0]['path_log_score'] <= done[-1]['path_log_score']:
            counts['score_bound_stopped_live'] = len(live)
            break
    results = []
    for rank, row in enumerate(_best(list(completed.values()), completed_limit)):
        if not math.isfinite(row['path_log_score']):
            raise ValueError('Nonfinite completed path score')
        results.append({'query_index': index, 'candidate_rank': rank,
                        'query_logit': float(window.query_logits[index]),
                        'chain': row['chain'], 'path_log_score': row['path_log_score'],
                        'terminated': True, 'stop_reason': 'canonical_stop'})
    return {'query_index': index, 'query_logit': float(window.query_logits[index]),
            'completed': results, 'counts': counts}


@torch.no_grad()
def decode_query_beam(model, features, sequence, max_exons, beam_width,
                      completed_limit, valid_bases=None):
    """Explicit-budget candidate search; not a substitute for frozen evaluation.

    Terminal paths and live histories have separate bounded pools. More than
    one returned chain per query requires a future chain-quality readout before
    any final annotation protocol, not duplicated existence gains by default.
    """
    if min(max_exons, beam_width, completed_limit) < 1:
        raise ValueError('Positive explicit exon/beam/completion limits required')
    window = model.for_window(features, valid_bases)
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Oriented DNA and valid feature lengths must match')
    if not bool(torch.isfinite(window.query_logits).all()):
        raise ValueError('Nonfinite query existence logit')
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    return [_beam_query(window, i, starts, acceptors, events, orf,
                        max_exons, beam_width, completed_limit)
            for i in range(model.num_queries)]
