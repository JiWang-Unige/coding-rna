"""A gene-query decoding with the existing canonical CDS grammar.

DNA and frozen features only; no B1 pool or reference injection. All legal
endpoint events remain available to each query. Greedy histories are neither
beam search nor a globally optimal chain decoder. This is an untrained oriented
window interface, not a complete FASTA annotation pipeline.
"""
import torch
from torch.nn import functional as F
from .candidates import Budget, ORFIndex, STOP
from .query_pointer import STAGES


def _end_events(sequence):
    return sorted([(i, 'donor') for i in range(1, len(sequence)-1)
                   if sequence[i:i+2] in ('GT', 'GC')] +
                  [(i+3, 'stop') for i in range(len(sequence)-2)
                   if sequence[i:i+3] in STOP])


def _legal_ends(orf, events, a, carry, first_exon):
    choices = []
    for b, kind in events:
        if b <= a or ((first_exon or kind == 'stop') and b-a < 3):
            continue
        remainder = orf.extend(a, b, carry, terminal=kind == 'stop')
        if remainder is not None:
            choices.append((b, kind, remainder))
    return choices


def _decode_query(window, index, starts, acceptors, events, orf, max_exons):
    model = window.model
    state, previous_end, carry, chain = window.initial[index], -1, '', []

    def result(terminated, reason):
        return {'query_index': index, 'query_logit': float(window.query_logits[index]),
                'chain': tuple(chain), 'terminated': terminated, 'stop_reason': reason}

    for j in range(max_exons):
        stage = 'start' if j == 0 else 'next_start'
        lower = 0 if j == 0 else previous_end + Budget().min_intron_bp
        entries = starts if j == 0 else [a for a in acceptors if a >= lower]
        if not entries:
            return result(False, 'no_start_motif' if j == 0 else 'no_acceptor')
        coordinates = torch.tensor(entries, device=state.device, dtype=torch.long)
        logits = window.pointer_logits(state, stage, lower)
        a = entries[int(logits[coordinates].argmax())]
        state = window.advance(state, stage, a)
        choices = _legal_ends(orf, events, a, carry, first_exon=j == 0)
        if not choices:
            return result(False, 'no_legal_exon_end')

        # Joint boundary/EOS choice on the DNA-legal events. EOS is conditioned
        # on each prospective post-boundary state, not on the preceding exon.
        coordinates = torch.tensor([b for b, _, _ in choices], device=state.device)
        tokens = (window.end_features[coordinates] + window.positions[coordinates] +
                  model.event.weight[STAGES['end']])
        states = model.recurrent(tokens, state.expand(len(choices), -1))
        eos = model.terminate(states).squeeze(-1)
        terminal = torch.tensor([kind == 'stop' for _, kind, _ in choices], device=state.device)
        branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
        boundary = F.log_softmax(window.pointer_logits(state, 'end', a), dim=-1)
        selected = int((boundary[coordinates] + branch).argmax())
        b, kind, carry = choices[selected]
        state = states[selected]
        chain.append((a, b))
        if kind == 'stop':
            return result(True, 'canonical_stop')
        previous_end = b
    # A cap never repairs a selected donor into a lower-scoring STOP.
    return result(False, 'exon_limit')


@torch.no_grad()
def decode_query_chains(model, features, sequence, max_exons, valid_bases=None):
    """Return raw per-query chains, including unterminated ones.

    Like B1: contiguous ATG start, GT/GC--AG, >=5bp introns, actual codon carry,
    no ambiguous coding bases/internal STOP, contiguous canonical final STOP.
    Split start/stop codons and other genetic codes are outside this grammar.
    Existence thresholding, duplicate/query matching and genome ownership are
    deliberately not chosen here; query logits are not calibrated probabilities.
    """
    if max_exons < 1:
        raise ValueError('A positive explicit exon limit is required')
    window = model.for_window(features, valid_bases)
    sequence = sequence.upper()
    if len(sequence) != window.length:
        raise ValueError('Oriented DNA and valid feature lengths must match')
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    events, orf = _end_events(sequence), ORFIndex(sequence)
    return [_decode_query(window, i, starts, acceptors, events, orf, max_exons)
            for i in range(model.num_queries)]
