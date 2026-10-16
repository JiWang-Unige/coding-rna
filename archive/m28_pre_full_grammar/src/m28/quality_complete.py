"""Paired complete-pool readouts on the SAME continued recurrent beam.

Live histories retain A5 cumulative scores. Complete paths are not discarded
before a frozen complete-chain scorer sees them. No raw-score early bound is
used: it cannot upper-bound future non-additive chain quality. This is not a
globally optimal search, prefix-quality model or reference-conditioned caller.
"""
import math
import torch
from torch.nn import functional as F
from .candidates import Budget, ORFIndex
from .query_beam import _best
from .query_decode import _end_events, _legal_ends
from .query_pointer import STAGES


def _continued_query(window, index, starts, acceptors, events, orf, max_exons, beam_width):
    live = [dict(state=window.initial[index], chain=(), carry='', previous_end=-1, path_log_score=0.)]
    completed = {}
    counts = dict(depths=0, entry_candidates=0, entry_pruned=0, legal_end_events=0,
                  continuation_pruned=0, live_at_exon_cap=0, completed_considered=0)
    for depth in range(max_exons):
        entries, entry_count = [], 0
        stage = 'start' if depth == 0 else 'next_start'
        for row in live:
            lower = 0 if depth == 0 else row['previous_end']+Budget().min_intron_bp
            coordinates = starts if depth == 0 else [a for a in acceptors if a >= lower]
            if not coordinates:
                continue
            logp = F.log_softmax(window.pointer_logits(row['state'], stage, lower), dim=-1)
            scores = logp[coordinates]
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal entry scores')
            order = torch.argsort(scores, descending=True, stable=True)[:beam_width]
            entry_count += len(coordinates)
            for k in order.tolist():
                a = coordinates[k]
                entries.append(dict(row, entry=a, state=window.advance(row['state'], stage, a),
                                    path_log_score=row['path_log_score']+float(scores[k])))
        kept = _best(entries, beam_width)
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
                      window.model.event.weight[STAGES['end']])
            states = window.model.recurrent(tokens, row['state'].expand(len(choices), -1))
            eos = window.model.terminate(states).squeeze(-1)
            terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices], device=eos.device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
            boundary = F.log_softmax(window.pointer_logits(row['state'], 'end', a), dim=-1)
            scores = boundary[coordinates]+branch
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal boundary/EOS scores')
            for k, (b, kind, carry) in enumerate(choices):
                chain, score = row['chain']+((a, b),), row['path_log_score']+float(scores[k])
                if kind == 'stop':
                    counts['completed_considered'] += 1
                    # Completed tensors are not retained; only legal geometry/score.
                    old = completed.get(chain)
                    if old is None or score > old['path_log_score']:
                        completed[chain] = dict(chain=chain, path_log_score=score)
                else:
                    continuations.append(dict(state=states[k], chain=chain, carry=carry,
                                              previous_end=b, path_log_score=score))
        live = _best(continuations, beam_width)
        counts['continuation_pruned'] += len(continuations)-len(live)
        counts['depths'] = depth+1
        if not live:
            break
        if depth == max_exons-1:
            counts['live_at_exon_cap'] = len(live)
    return dict(query_index=index, query_logit=float(window.query_logits[index]),
                considered=list(completed.values()), counts=counts)


@torch.no_grad()
def decode_query_complete_pool(model, features, sequence, score_completed,
                               max_exons, beam_width, completed_limit, valid_bases=None):
    """Score DNA-only complete chains once per geometry, then pair top-k readouts.

    Both pools share the same examined paths. Quality changes only completed
    retention; it does not change live prefix states/ranking. Nonpositive quality
    chains remain free candidates. Existence is metadata, never a chain gain.
    """
    if min(max_exons, beam_width, completed_limit) < 1:
        raise ValueError('Positive explicit exon/beam/completion limits required')
    window = model.for_window(features, valid_bases)
    sequence = sequence.upper()
    if len(sequence) != window.length or not bool(torch.isfinite(window.query_logits).all()):
        raise ValueError('Aligned DNA/features and finite query logits required')
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    queries = [_continued_query(window, q, starts, acceptors, events, orf, max_exons, beam_width)
               for q in range(model.num_queries)]
    chains = sorted({r['chain'] for q in queries for r in q['considered']})
    qualities = list(score_completed(chains)) if chains else []
    if len(qualities) != len(chains) or not all(math.isfinite(float(v)) for v in qualities):
        raise ValueError('One finite frozen quality score per complete geometry required')
    scores = {chain:float(score) for chain, score in zip(chains, qualities)}
    for query in queries:
        candidates = [dict(r, quality_logit=scores[r['chain']], query_index=query['query_index'],
                           terminated=True, stop_reason='canonical_stop') for r in query.pop('considered')]
        query['path_completed'] = _best(candidates, completed_limit)
        query['quality_completed'] = sorted(candidates, key=lambda r:(-r['quality_logit'], r['chain']))[:completed_limit]
    return dict(queries=queries, all_completed=[dict(chain=c, quality_logit=scores[c]) for c in chains],
                early_score_bound_used=False, quality_used_for_live_history=False)
