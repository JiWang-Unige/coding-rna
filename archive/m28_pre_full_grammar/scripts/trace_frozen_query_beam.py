"""Observe frozen A5 beam decisions without injecting or retaining reference paths.

This one-off diagnostic mirrors query_beam._beam_query. Its returned native
result must be reconciled against that original implementation on actual data.
Reference geometry is used only after each ordinary pool has been constructed.
"""
import math
import torch
from torch.nn import functional as F
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.query_pointer import STAGES
from src.m28.query_beam import _best, _beam_query


def view(row):
    return dict(chain=row['chain'], entry=row.get('entry'),
                path_log_score=row['path_log_score'])


def pool_loss(stage, depth, target, pool, kept):
    ordered = _best(pool, len(pool))
    return dict(stage=stage, exon=depth+1, target=view(target),
                target_rank=next(i+1 for i, r in enumerate(ordered) if r is target),
                pool_size=len(pool), retained_size=len(kept),
                cutoff=view(kept[-1]) if kept else None)


@torch.inference_mode()
def trace_query(window, sequence, index, references, max_exons=128,
                beam_width=4, completed_limit=4):
    """Return unchanged natural completions and each target's first actual loss.

No forced prefix is advanced. A reference never enters a candidate pool and
cannot change a ranking, normalizer, budget, EOS decision or stopping rule.
"""
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Original oriented DNA and valid feature lengths required')
    targets = {key: tuple(tuple(x) for x in chain) for key, chain in references.items()}
    if any(not chain or len(chain) > max_exons for chain in targets.values()):
        raise ValueError('Only supported nonempty chains within the frozen cap')
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    model = window.model
    counts = dict(entry_candidates=0, entry_pruned=0, legal_end_events=0,
                  completed_considered=0, continuation_pruned=0,
                  live_at_exon_cap=0, score_bound_stopped_live=0)
    losses = {key: None for key in targets}
    live = [dict(state=window.initial[index], chain=(), carry='',
                 previous_end=-1, path_log_score=0.)]
    completed = {}
    for depth in range(max_exons):
        entries = []; depth_candidates = 0
        for row in live:
            lower = 0 if depth == 0 else row['previous_end'] + Budget().min_intron_bp
            coordinates = starts if depth == 0 else [a for a in acceptors if a >= lower]
            if not coordinates:
                continue
            stage = 'start' if depth == 0 else 'next_start'
            logprob = F.log_softmax(window.pointer_logits(row['state'], stage, lower), dim=-1)
            coordinate_tensor = torch.tensor(coordinates, device=logprob.device)
            scores = logprob[coordinate_tensor]
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal entry score')
            order = torch.argsort(scores, descending=True, stable=True)[:beam_width]
            counts['entry_candidates'] += len(coordinates)
            depth_candidates += len(coordinates)
            kept_coordinates = [coordinates[k] for k in order.tolist()]
            for key, chain in targets.items():
                if losses[key] is None and depth < len(chain) and row['chain'] == chain[:depth]:
                    a = chain[depth][0]
                    if a not in coordinates:
                        raise ValueError('A grammar-supported target has no legal entry')
                    if a not in kept_coordinates:
                        k = coordinates.index(a)
                        values = scores.tolist()
                        rank = 1 + sum(v > values[k] or (v == values[k] and x < a)
                                       for x, v in zip(coordinates, values))
                        losses[key] = dict(stage='local_'+stage+'_topk', exon=depth+1,
                            target=dict(chain=row['chain'], entry=a,
                                        path_log_score=row['path_log_score']+float(scores[k])),
                            target_rank=rank, pool_size=len(coordinates),
                            retained_size=len(order), cutoff_coordinate=kept_coordinates[-1])
            for k in order.tolist():
                a = coordinates[k]
                entries.append({**row, 'entry': a,
                                'state': window.advance(row['state'], stage, a),
                                'path_log_score': row['path_log_score'] + float(scores[k])})
        kept_entries = _best(entries, beam_width)
        counts['entry_pruned'] += depth_candidates - len(kept_entries)
        for key, chain in targets.items():
            if losses[key] is None and depth < len(chain):
                target = next((r for r in entries if r['chain'] == chain[:depth]
                               and r['entry'] == chain[depth][0]), None)
                if target is not None and not any(r is target for r in kept_entries):
                    losses[key] = pool_loss('global_entry_pool_pruned', depth, target, entries, kept_entries)
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
            terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices], device=row['state'].device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
            boundary = F.log_softmax(window.pointer_logits(row['state'], 'end', a), dim=-1)
            scores = boundary[coordinates] + branch
            if not bool(torch.isfinite(scores).all()):
                raise ValueError('Nonfinite legal end/EOS score')
            for k, (b, kind, carry) in enumerate(choices):
                new = dict(state=states[k], chain=row['chain']+((a, b),), carry=carry,
                           previous_end=b, path_log_score=row['path_log_score']+float(scores[k]))
                if kind == 'stop':
                    counts['completed_considered'] += 1
                    old = completed.get(new['chain'])
                    if old is None or new['path_log_score'] > old['path_log_score']:
                        completed[new['chain']] = new
                else:
                    continuations.append(new)
        full_done = list(completed.values())
        done = _best(full_done, completed_limit)
        completed = {r['chain']: r for r in done}
        live = _best(continuations, beam_width)
        counts['continuation_pruned'] += len(continuations) - len(live)
        for key, chain in targets.items():
            if losses[key] is not None:
                continue
            if depth+1 < len(chain):
                target = next((r for r in continuations if r['chain'] == chain[:depth+1]), None)
                if target is not None and not any(r is target for r in live):
                    losses[key] = pool_loss('continuation_pool_pruned', depth, target, continuations, live)
            else:
                target = next((r for r in full_done if r['chain'] == chain), None)
                if target is not None and chain not in completed:
                    losses[key] = pool_loss('completion_pool_pruned', depth, target, full_done, done)
        reason = None
        if not live:
            reason = 'no_live_histories'
        elif depth == max_exons-1:
            counts['live_at_exon_cap'] = len(live); reason = 'exon_cap'
        elif len(done) == completed_limit and live[0]['path_log_score'] <= done[-1]['path_log_score']:
            counts['score_bound_stopped_live'] = len(live); reason = 'retained_beam_score_bound'
        if reason is not None:
            for key, chain in targets.items():
                if losses[key] is None and chain not in completed:
                    target = next((r for r in live if r['chain'] == chain[:depth+1]), None)
                    if target is None:
                        raise ValueError('Supported target disappeared without an observed pruning event')
                    losses[key] = dict(stage=reason, exon=depth+1, target=view(target),
                                       completion_cutoff=view(done[-1]) if done else None)
            break
    results = []
    for rank, row in enumerate(_best(list(completed.values()), completed_limit)):
        if not math.isfinite(row['path_log_score']):
            raise ValueError('Nonfinite completed path score')
        results.append(dict(query_index=index, candidate_rank=rank,
            query_logit=float(window.query_logits[index]), chain=row['chain'],
            path_log_score=row['path_log_score'], terminated=True, stop_reason='canonical_stop'))
    natural = dict(query_index=index, query_logit=float(window.query_logits[index]), completed=results, counts=counts)
    traces = {key: dict(returned=chain in completed, first_loss=losses[key]) for key, chain in targets.items()}
    if any(t['returned'] == (t['first_loss'] is not None) for t in traces.values()):
        raise ValueError('Returned chains and first-loss events disagree')
    return natural, traces


@torch.inference_mode()
def original_query(window, sequence, index, max_exons=128, beam_width=4, completed_limit=4):
    sequence = sequence.upper()
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    return _beam_query(window, index, starts, acceptors, _end_events(sequence), ORFIndex(sequence),
                       max_exons, beam_width, completed_limit)


def reconcile(original, replay):
    if original['counts'] != replay['counts'] or original['query_index'] != replay['query_index']:
        raise ValueError('Observed replay changed original query identity or search counts')
    if original['query_logit'] != replay['query_logit'] or len(original['completed']) != len(replay['completed']):
        raise ValueError('Observed replay changed original query existence or completion count')
    for a, b in zip(original['completed'], replay['completed']):
        if ({k:v for k,v in a.items() if k != 'path_log_score'} !=
                {k:v for k,v in b.items() if k != 'path_log_score'} or
                not math.isclose(a['path_log_score'], b['path_log_score'], rel_tol=0., abs_tol=1e-3)):
            raise ValueError('Observed replay changed completion order, geometry or path score')
