"""Gradient-bearing replay of an already selected original-search path.

No candidate collection, rollout, new decoder or model is implemented here.
The caller must collect current search decisions and independently qualify
cutoffs; an eligible primary or unresolved cutoff is never a negative.
"""
import math
import torch
from torch.nn import functional as F
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.query_prefix_value import START, NEXT_START, DONOR, STOP


def replay_node_score(window, sequence, query_index, node):
    """Recompute the native cumulative geometric/EOS score with a live graph.

All closed exons of an entry or DONOR node use the nonterminal EOS branch;
only the last exon of a STOP node uses the terminal branch. Intron lower
bounds and full geometric normalizers are the original raw-search ones,
not the teacher's previous_end+1 or a motif-only normalization.
"""
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Original oriented DNA and valid window must agree')
    chain = tuple(tuple(e) for e in node['chain'])
    entry, event = node['entry'], node['event']
    if event in (START, NEXT_START):
        if entry is None or event != (NEXT_START if chain else START):
            raise ValueError('Entry event must match the observed history')
    elif event in (DONOR, STOP):
        if not chain or entry is not None:
            raise ValueError('Closed exon event requires a nonempty chain and no pending entry')
    else:
        raise ValueError('Only original entry, donor and stop nodes are supported')
    if chain:
        window.validate_chain(chain)
    orf, ends = ORFIndex(sequence), _end_events(sequence)
    state, carry, previous_end = window.initial[query_index], '', -1
    score = state.sum()*0

    def enter(state, a, depth, previous_end):
        stage = 'start' if depth == 0 else 'next_start'
        lower = 0 if depth == 0 else previous_end+Budget().min_intron_bp
        legal = (0 <= a < window.length and a >= lower and
                 (sequence[a:a+3] == 'ATG' if depth == 0 else sequence[a-2:a] == 'AG'))
        if not legal:
            raise ValueError('Selected entry is outside the original DNA search grammar')
        logp = F.log_softmax(window.pointer_logits(state, stage, lower), dim=-1)[a]
        return window.advance(state, stage, a), logp

    for depth, (a, b) in enumerate(chain):
        state, logp = enter(state, a, depth, previous_end)
        kind = 'stop' if event == STOP and depth == len(chain)-1 else 'donor'
        choices = _legal_ends(orf, ends, a, carry, first_exon=depth == 0)
        selected = [(x, k, c) for x, k, c in choices if (x, k) == (b, kind)]
        if len(selected) != 1:
            raise ValueError('Selected end/EOS branch is outside the original ORF grammar')
        boundary = F.log_softmax(window.pointer_logits(state, 'end', a), dim=-1)[b]
        state = window.advance(state, 'end', b)
        eos = window.model.terminate(state).squeeze(-1)
        branch = F.logsigmoid(eos if kind == 'stop' else -eos)
        score = score+logp+boundary+branch
        carry, previous_end = selected[0][2], b
    if entry is not None:
        _, logp = enter(state, entry, len(chain), previous_end)
        score = score+logp
    return score


def qualified_retention_margin(target_score, cutoff_score, cutoff_label, margin):
    """A score hinge only for a separately qualified known-nonmatch cutoff.

Margin is explicit; no formal training coefficient or sampling protocol is
chosen here. A zero loss on positive/protected cutoffs does not mean no other
safe competitor exists, since a single cutoff is not the whole search pool.
"""
    if cutoff_label not in (-1, 0, 1) or not math.isfinite(margin) or margin < 0:
        raise ValueError('Explicit 1/0/-1 qualification and finite nonnegative margin required')
    if target_score.ndim != 0 or cutoff_score.ndim != 0:
        raise ValueError('One same-query target and cutoff score required')
    if cutoff_label != 0:
        return (target_score+cutoff_score)*0
    return F.relu(margin+cutoff_score-target_score)
