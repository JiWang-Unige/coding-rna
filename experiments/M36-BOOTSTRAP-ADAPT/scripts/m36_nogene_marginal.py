"""Per-base 'no gene covers this base' marginal under the frozen complete-gene-set grammar (protocol 22 §2 S1 pseudo-background; §7 item 4).

In the frozen automaton (`full_grammar_boundary_messages._transitions`) the background self-loop B -> B is the only transition that emits role 0,
and every gene occupies a B -> ... -> B excursion.  Hence, with alpha_B(i) = log mass of all prefixes of length i that end in the background state and
beta_B(j) = log mass of all suffixes starting at j in the background state (the quantities the frozen `_prefix_queries` / `_suffix_queries` already
report at query positions),      log P(base i not covered by any gene) = alpha_B(i) + gain[i, 0] + beta_B(i + 1) - log Z,   gain[i, 0] = 0.
This module only adds the query at every boundary 0..n (one forward + one backward pass, the same cost as the M33 gene posterior) using the frozen
scaled-mass kernels; no frozen file is edited, no gene set is pruned, no threshold enters the computation.
usage:  m36_nogene_marginal.py run <run_dir> <seq> <strand> <tau> <out_dir> [--fasta genome.fa]
            reads <run_dir>/pot/<seq>_<strand>.npy (oriented [L,15] float32, as m33_infer writes) ; writes <out_dir>/tau<tau>/nogene_<seq>_<strand>.npy
            (float32 P(no gene), oriented) and <out_dir>/tau<tau>/nogene_<seq>_<strand>.json (log Z, pseudo-background segments at P >= 0.99, seconds)
        m36_nogene_marginal.py bench <length_bp> [seed]      synthetic timing (random DNA + random potentials), also checks log Z against the frozen posterior
"""
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
import torch

import m36_common as W
from src.m28.full_grammar_boundary_compiled import ALPHABET, _state_table
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries


def no_gene_marginal(sequence, logits, max_exons=W.MAX_EXONS, block_bp=256):
    """sequence: oriented DNA (str); logits: [n,15] CPU float32/64 energies (tau already applied).
    -> dict(log_partition, log_p_nogene [n] float64, p_nogene [n] float64, alpha_B [n+1], beta_B [n+1])"""
    sequence = sequence.upper()
    n = len(sequence)
    if (not n or tuple(logits.shape) != (n, 15) or logits.device.type != 'cpu' or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32/FP64 roles aligned to nonempty oriented DNA required')
    states, destinations, roles, background = _state_table(max_exons)
    codes = np.array([ALPHABET.find(b) if b in ALPHABET else 4 for b in sequence], dtype=np.int8)
    gain = (logits.detach().double() - logits.detach()[:, :1].double()).numpy()
    positions = np.arange(n + 1, dtype=np.int64)
    log_z, alpha = _scaled_prefix_queries(codes, gain, positions, destinations, roles, background, block_bp)
    reverse_z, beta = _scaled_suffix_queries(codes, gain, positions, destinations, roles, background, block_bp)
    if not math.isclose(log_z, reverse_z, rel_tol=1e-9, abs_tol=1e-6) or log_z < -1e-8:
        raise ValueError(f'Forward/backward total mass differs: {log_z} vs {reverse_z}')
    logp = alpha[:-1] + beta[1:] - log_z                      # + gain[:, 0] which is identically 0
    if np.nanmax(logp) > 1e-6:
        raise ValueError('event mass exceeds the partition')
    logp = np.minimum(logp, 0.)
    return dict(log_partition=float(log_z), reverse_log_partition=float(reverse_z), log_p_nogene=logp, p_nogene=np.exp(logp), alpha_B=alpha, beta_B=beta,
                state_count=len(states), max_exons=max_exons, reference_used=False, posterior_threshold_applied=False)


def segments(p, threshold=W.PSEUDO_BACKGROUND_MIN_P, min_len=1):
    """oriented runs with p >= threshold -> [(x, y)] (0-based half-open, oriented frame)"""
    hit = np.asarray(p) >= threshold
    d = np.diff(np.concatenate([[0], hit.astype(np.int8), [0]]))
    return [(int(x), int(y)) for x, y in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)) if y - x >= min_len]


def to_genomic(intervals, strand, length):
    return list(intervals) if strand == '+' else [(length - y, length - x) for x, y in reversed(intervals)]


def run(argv):
    run_dir, seq, strand, tau, outd = Path(argv[0]), argv[1], argv[2], float(argv[3]), Path(argv[4])
    fasta = argv[argv.index('--fasta') + 1] if '--fasta' in argv else None
    torch.set_num_threads(1)
    out = outd / f'tau{tau:g}' / f'nogene_{seq}_{strand}'
    if Path(str(out) + '.json').exists():
        print('exists', out)
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    pot = np.load(run_dir / 'pot' / f'{seq}_{strand}.npy', mmap_mode='r')
    L = pot.shape[0]
    if fasta:
        import manifest as MF
        dna = MF.read_allowed(Path(fasta), {seq: L})[seq]
    else:
        M29 = W.M30.register()
        dna = M29.read_sequence(seq)
    import manifest as MF
    dna = re.sub('[^ACGT]', 'N', dna.upper())
    assert len(dna) == L, (len(dna), L)
    dna = dna if strand == '+' else MF.rc(dna)
    t0 = time.perf_counter()
    lg = torch.from_numpy(np.ascontiguousarray(pot * np.float32(tau)))
    r = no_gene_marginal(dna, lg)
    p = r['p_nogene'].astype(np.float32)
    np.save(str(out) + '.npy.tmp.npy', p)
    Path(str(out) + '.npy.tmp.npy').replace(str(out) + '.npy')
    seg = segments(p)
    json.dump(dict(seq=seq, strand=strand, tau=tau, length=L, log_partition=r['log_partition'], threshold=W.PSEUDO_BACKGROUND_MIN_P,
                   pseudo_background_oriented=seg, pseudo_background_genomic=to_genomic(seg, strand, L),
                   covered_bp=int(sum(y - x for x, y in seg)), mean_p_nogene=float(p.mean()), seconds=time.perf_counter() - t0),
              open(str(out) + '.json.tmp', 'w'))
    Path(str(out) + '.json.tmp').replace(str(out) + '.json')
    print(seq, strand, f'tau={tau:g}', L, f'{time.perf_counter() - t0:.0f}s', 'pseudo-bg bp', sum(y - x for x, y in seg), flush=True)


def synthetic(length, seed=0, bg=0.5):
    """random DNA with a few planted ATG...GT..AG...stop structures and random N(0,1) potentials with a background bias"""
    rng = np.random.default_rng(seed)
    dna = rng.choice(list('ACGT'), length)
    for k in range(max(1, length // 3000)):
        i = int(rng.integers(0, max(1, length - 200)))
        gene = 'ATG' + ''.join(rng.choice(['GCT', 'GAA', 'AAA', 'CCT', 'GGC', 'TTC'], 12)) + 'GTAAGTTTTCAG' + ''.join(rng.choice(['GCT', 'GAA', 'AAA'], 9)) + 'TAA'
        dna[i:i + len(gene)] = list(gene[:length - i])
    dna = ''.join(dna)
    lg = rng.normal(0, 1, (length, 15)).astype(np.float32)
    lg[:, 0] += bg
    return dna, torch.from_numpy(lg)


def bench(argv):
    length = int(float(argv[0]))
    seed = int(argv[1]) if len(argv) > 1 else 0
    torch.set_num_threads(1)
    dna, lg = synthetic(min(length, 20000), seed)
    no_gene_marginal(dna, lg)                                           # numba compile warm-up
    dna, lg = synthetic(length, seed)
    t0 = time.perf_counter()
    r = no_gene_marginal(dna, lg)
    dt = time.perf_counter() - t0
    out = dict(length=length, seconds=round(dt, 2), seconds_per_Mb=round(dt / (length / 1e6), 2), projected_20Mb_minutes=round(dt / (length / 1e6) * 20 / 60, 1),
               state_count=r['state_count'], log_partition=r['log_partition'], mean_p_nogene=float(r['p_nogene'].mean()),
               pseudo_bg_fraction_at_p99=float(np.mean(r['p_nogene'] >= 0.99)))
    if length <= 300000:
        import src.m28.full_grammar_gene_posterior as G
        t1 = time.perf_counter()
        z = G.complete_gene_set_posteriors(dna, lg, [], W.MAX_EXONS)['log_partition']
        out.update(frozen_log_partition=z, frozen_seconds=round(time.perf_counter() - t1, 2), log_partition_abs_diff=abs(z - r['log_partition']))
    print(json.dumps(out), flush=True)
    return out


if __name__ == '__main__':
    W.K.repoint()
    if sys.argv[1] == 'run':
        run(sys.argv[2:])
    elif sys.argv[1] == 'bench':
        bench(sys.argv[2:])
