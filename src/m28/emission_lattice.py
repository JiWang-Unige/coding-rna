"""Source-only exact oriented gene-set decoder for the existing CDS grammar.

Optimizes the additive 15-role-minus-intergenic objective, not a query model,
calibrated likelihood, or unrestricted biological grammar. All motif-legal
edges are considered. Floating-point accumulation is FP64, not exact arithmetic.
"""
from collections import defaultdict
from dataclasses import dataclass

import torch

from .candidates import Budget, ORFIndex, STOP


class _Margins:
    def __init__(self, logits):
        gain = (logits.double() - logits[:, :1].double()).cpu()
        self.gain = gain.tolist()
        self.prefix = torch.cat((torch.zeros(1, 15, dtype=torch.float64),
                                 gain.cumsum(0))).T.tolist()
        positions = torch.arange(len(gain))
        coding = torch.stack([gain[positions, 1 + (positions + offset) % 3]
                              for offset in range(3)], dim=1)
        self.coding = torch.cat((torch.zeros(1, 3, dtype=torch.float64),
                                 coding.cumsum(0))).T.tolist()

    def interval(self, a, b, role):
        return self.prefix[role][b] - self.prefix[role][a]

    def coding_interval(self, a, b, phase):
        offset = (phase - a) % 3
        return self.coding[offset][b] - self.coding[offset][a]

    def exon_gain(self, a, b, phase, first, terminal):
        """Includes donor now; its future intron is charged at the acceptor."""
        size = b - a
        value = self.coding_interval(a, b, phase)
        if first:
            value += self.interval(a, a + 3, 13) - self.coding_interval(a, a + 3, 0)
        elif not (terminal and size == 3):
            value += self.gain[a][10 + phase] - self.gain[a][1 + phase]
        if terminal:
            value += (self.interval(b - 3, b, 14) -
                      self.coding_interval(b - 3, b, (phase + size - 3) % 3))
            if first and size == 3:
                value -= self.interval(a, b, 13) - self.coding_interval(a, b, 0)
        elif not (first and size == 3):
            end_phase = (phase + size - 1) % 3
            prior = 10 + phase if not first and size == 1 else 1 + end_phase
            value += self.gain[b - 1][7 + end_phase] - self.gain[b - 1][prior]
        return value


@dataclass(frozen=True)
class _Solution:
    score: float = 0.
    chains: tuple = ()
    gains: tuple = ()


@dataclass(frozen=True)
class _Live:
    score: float
    preceding: _Solution
    chain: tuple
    carry: str


def _retain(table, key, row):
    previous = table.get(key)
    if previous is None or row.score > previous.score:
        table[key] = row
        return previous is not None
    return True


@torch.no_grad()
def decode_emission_lattice(sequence, logits, max_exons):
    """Best nonoverlapping complete-CDS set in one oriented, valid-length window.

Raw support matches the query grammar: contiguous ATG and final STOP,
GT/GC--AG introns of at least five bases, actual unfinished codons, no ambiguous
coding bases or internal STOP, and the explicit exon cap. Intronic ambiguity
is allowed. No site top-k, query/pointer, beam, labels, or candidate pool enters.

The state contains exon count and actual codon carry. Mature donors are merged
only after subtracting the same intron-role prefix, making their future scores
identical up to a constant. Terminal chains update the intergenic optimum;
adjacent genes are legal. Equal-scoring completions do not replace that optimum,
so an all-zero objective returns null. Other ties retain deterministic first
arrival rather than promising the original finite-pool DP's identity policy.
"""
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if logits.ndim != 2 or tuple(logits.shape) != (length, 15) or not length:
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if not bool(torch.isfinite(logits).all()):
        raise ValueError('Finite local potentials required')
    scores, orf = _Margins(logits), ORFIndex(sequence)
    starts = {i for i in range(length - 2) if sequence[i:i + 3] == 'ATG'}
    acceptors = {i + 2 for i in range(length - 2) if sequence[i:i + 2] == 'AG'}
    donors = {i for i in range(1, length - 1) if sequence[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(length - 2) if sequence[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    maturities = {d + gap: d for d in donors if d + gap <= length}
    positions = sorted(starts | acceptors | stops | set(maturities))
    pending, active, completed, edge_cache = defaultdict(dict), {}, {}, {}
    best = _Solution()
    counts = dict(legal_exon_edges=0, exon_extensions=0, donor_states_created=0,
                  donor_states_dominated=0, intron_states_dominated=0,
                  completed_extensions=0, exon_cap_nonterminal_edges=0,
                  maximum_active_states=0)

    def legal_ends(a, carry, first):
        key = (a, carry, first)
        if key not in edge_cache:
            choices = []
            for b, terminal in events:
                if b <= a or ((first or terminal) and b - a < 3):
                    continue
                remainder = orf.extend(a, b, carry, terminal)
                if remainder is not None:
                    # The increment depends on entry/carry/first, not the
                    # preceding gene set or this history's exon count.
                    delta = scores.exon_gain(a, b, len(carry), first, terminal)
                    choices.append((b, terminal, remainder, delta))
            edge_cache[key] = choices
            counts['legal_exon_edges'] += len(choices)
        return edge_cache[key]

    for a in positions:
        solution = completed.pop(a, None)
        if solution is not None and solution.score > best.score:
            best = solution
        d = maturities.get(a)
        if d is not None:
            for key, row in sorted(pending.pop(d, {}).items()):
                intron_role = 4 + (len(row.carry) - 1) % 3
                normalized = _Live(row.score - scores.prefix[intron_role][d],
                                   row.preceding, row.chain, row.carry)
                counts['intron_states_dominated'] += _retain(active, key, normalized)
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        entries = []
        if a in starts:
            entries.append(_Live(best.score, best, (), ''))
        if a in acceptors:
            for _, row in sorted(active.items()):
                intron_role = 4 + (len(row.carry) - 1) % 3
                entries.append(_Live(row.score + scores.prefix[intron_role][a],
                                     row.preceding, row.chain, row.carry))
        for row in entries:
            first = not row.chain
            for b, terminal, carry, delta in legal_ends(a, row.carry, first):
                counts['exon_extensions'] += 1
                score = row.score + delta
                exon_count = len(row.chain) + 1
                if terminal:
                    counts['completed_extensions'] += 1
                    old = completed.get(b)
                    if old is None or score > old.score:
                        chain = row.chain + ((a, b),)
                        completed[b] = _Solution(score, row.preceding.chains + (chain,),
                            row.preceding.gains + (score - row.preceding.score,))
                elif exon_count < max_exons:
                    counts['donor_states_created'] += 1
                    table, key = pending[b], (exon_count, carry)
                    old = table.get(key)
                    counts['donor_states_dominated'] += old is not None
                    if old is None or score > old.score:
                        # Materialize only retained histories, preserving the
                        # same strict comparison and first-arrival ties.
                        table[key] = _Live(score, row.preceding, row.chain + ((a, b),), carry)
                else:
                    counts['exon_cap_nonterminal_edges'] += 1
    return dict(chains=best.chains, gain=best.score, chain_gains=best.gains,
                counts=counts, reference_used=False)
