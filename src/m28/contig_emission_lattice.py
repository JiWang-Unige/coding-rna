"""Joint contig grammar decoding under a NEW base-owned context objective.

Each base takes all fifteen FP32 roles from one predetermined containing
window. This is not equivalent to old whole-chain owner-window scoring.
No model, FASTA file, label, candidate pool, or fitted correction is loaded.
"""
import torch

from .c0_decode import genomic_chain
from .emission_lattice import decode_emission_lattice


def context_segments(length, windows):
    """Genomic half-open segments owned by nearest containing window center.

Distance is measured from the base center p + 1/2. Equal distances choose the
lower genomic window start on BOTH strands. Supported windows cover the contig
and have strictly increasing starts and ends, as does the original M28 grid.
"""
    windows = [tuple(w) for w in windows]
    if (not length or not windows or windows[0][0] != 0 or windows[-1][1] != length or
            any(not 0 <= a < b <= length for a, b in windows)):
        raise ValueError('Nonempty complete contig window coverage required')
    cuts = [0]
    for (a, b), (c, d) in zip(windows, windows[1:]):
        if not (a < c <= b < d):
            raise ValueError('Increasing, nonnested windows without gaps required')
        # 4*p + 2 <= (a+b) + (c+d) belongs to the earlier window.
        # Clamp to the actual overlap for unequal valid window lengths.
        cuts.append(min(b, max(c, (a + b + c + d + 2) // 4)))
    cuts.append(length)
    return tuple((i, cuts[i], cuts[i + 1]) for i in range(len(windows)))


@torch.no_grad()
def stitch_contig_emissions(sequence, windows, strand, records):
    """Copy valid-length oriented window rows; never average or rescore them.

Records arrive once in increasing GENOMIC window order and contain
window_index, sequence and local_potentials. Their DNA must equal the relevant
oriented slice of the already normalized ACGTN contig. Padding is not accepted.
Minus-strand rows are already in local reverse-complement order: coordinates
change, but role columns and coding phases must NOT be reversed or permuted.
"""
    sequence = sequence.upper()
    if strand not in ('+', '-') or set(sequence) - set('ACGTN'):
        raise ValueError('One strand and normalized ACGTN DNA required')
    windows = [tuple(w) for w in windows]
    segments = context_segments(len(sequence), windows)
    rc = lambda dna: dna.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
    oriented = sequence if strand == '+' else rc(sequence)
    potentials = torch.empty(len(sequence), 15, dtype=torch.float32)
    received = 0
    for i, row in enumerate(records):
        if i >= len(windows) or row['window_index'] != i:
            raise ValueError('Exactly one record per window in genomic order required')
        a, b = windows[i]
        expected = sequence[a:b] if strand == '+' else rc(sequence[a:b])
        logits = row['local_potentials']
        if row['sequence'].upper() != expected:
            raise ValueError('Oriented window DNA does not match this contig slice')
        if (logits.device.type != 'cpu' or logits.dtype != torch.float32 or
                tuple(logits.shape) != (b - a, 15) or not bool(torch.isfinite(logits).all())):
            raise ValueError('Finite CPU FP32 valid-length fifteen-role rows required')
        _, left, right = segments[i]
        if strand == '+':
            potentials[left:right] = logits[left - a:right - a]
        else:
            potentials[len(sequence) - right:len(sequence) - left] = logits[b - right:b - left]
        received += 1
    if received != len(windows):
        raise ValueError('Missing contig window record')
    return dict(sequence=oriented, local_potentials=potentials,
                genomic_context_segments=segments)


@torch.no_grad()
def decode_contig_emission_lattice(sequence, windows, strand, records, max_exons):
    """One joint oriented-contig solve, not local selection then owner filtering.

Grammar, FP64 accumulation, role weights, null and ties are unchanged from the
window decoder. Context/score and allowed spans are different. Genomic chains
stay paired with oriented chain_gains; their ORDER on minus is transcriptional,
not increasing genomic start. Chromosome-scale cost remains unverified.
"""
    joined = stitch_contig_emissions(sequence, windows, strand, records)
    result = decode_emission_lattice(joined['sequence'], joined['local_potentials'], max_exons)
    result['genomic_chains'] = tuple(tuple(genomic_chain(chain, 0, len(sequence), strand))
                                     for chain in result['chains'])
    result['genomic_context_segments'] = joined['genomic_context_segments']
    result['context_policy'] = 'nearest_containing_genomic_window_center_per_base_lower_start_ties'
    return result
