"""CPU tests of the M41 intron-length grammar (no GPU, no data).
  * states OFF  : state table, logZ, gradient, coarse-channel loss and skyline decode equal the frozen M28/M29 code (bit-identical);
  * states ON   : logZ equals an explicit enumeration of every complete gene set with closed-form intron-length terms; the decoder's optimum equals the
                  enumerated maximum (max over latent paths), the decoded chains are legal; equal-rate branches reproduce the frozen logZ exactly;
  * gradients   : torch.autograd.gradcheck wrt the logits and wrt the three prior scalars.
Run: python -m pytest experiments/M41-INTRON-LENGTH/tests -q
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (str(HERE.parent / 'scripts'), str(ROOT / 'pipeline/m28'), str(ROOT / 'experiments/M29-CARBON-A-GAP/scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)

import m41_grammar as GR                                                     # noqa: E402
import m41_prior as P                                                        # noqa: E402
from src.m28.full_grammar_boundary_compiled import _state_table, compiled_boundary_full_grammar_log_partition   # noqa: E402
from src.m28.full_grammar_coarse_channel import coarse_channel_nll           # noqa: E402
from src.m28.full_grammar_partition import _TensorMargins                    # noqa: E402
from src.m28.candidates import Budget, ORFIndex, STOP                         # noqa: E402

ON1, ON3 = GR.IntronLength(True, 1), GR.IntronLength(True, 3)
TH = torch.tensor([1.2, 2.0, 4.0], dtype=torch.float64)                       # pi_S .77, short mean ~12 bp interior, long mean ~60: lengths that matter in short tests
P0 = 0.9


def rand_case(n, seed, scale=2.):
    rng = np.random.default_rng(seed)
    seq = ''.join(rng.choice(list('ACGT'), n))
    lg = torch.from_numpy((rng.standard_normal((n, 15)) * scale).astype(np.float32))
    return seq, lg


def planted(seed, n_units=3):
    """Short sequence rich in ATG / GT / AG / stops so that multi-exon genes with several intron lengths exist."""
    rng = np.random.default_rng(seed)
    parts = ['C', 'ATG']
    for _ in range(n_units):
        parts += [''.join(rng.choice(list('ACG'), int(rng.integers(1, 4)))), 'GT', ''.join(rng.choice(list('ACT'), int(rng.integers(1, 9)))), 'AG']
    parts += [''.join(rng.choice(list('CG'), int(rng.integers(0, 3)))), 'TAA', 'C', 'ATGCC', 'GTAAG', 'CTAG', 'GC']
    seq = ''.join(parts)
    lg = torch.from_numpy((rng.standard_normal((len(seq), 15)) * 1.5).astype(np.float64))
    return seq, lg


# ---- states OFF == frozen ----------------------------------------------------------------------------------------------------------
def test_off_state_table_is_frozen():
    for cap in (1, 3, 128):
        s0, d0, r0, b0 = _state_table(cap)
        s1, d1, r1, t1, b1 = GR.state_table(cap, GR.OFF.key())
        assert s0 == s1 and b0 == b1 and np.array_equal(d0, d1) and np.array_equal(r0, r1) and not t1.any()


def test_off_logz_and_gradient_bit_identical():
    seq, lg = rand_case(1300, 1)
    allowed = torch.ones((1300, 15), dtype=torch.bool)
    allowed[400:420] = False
    allowed[400:420, 0] = True
    for kw in (dict(), dict(allowed=allowed, cuts=(300, 900))):
        a = lg.clone().requires_grad_(True)
        r0 = compiled_boundary_full_grammar_log_partition(seq, a, 128, checkpoint_bp=256, **kw)
        r0['log_partition'].backward()
        b = lg.clone().requires_grad_(True)
        r1 = GR.log_partition(seq, b, 128, checkpoint_bp=256, **kw)
        r1['log_partition'].backward()
        assert r0['log_partition'].item() == r1['log_partition'].item()
        assert torch.equal(a.grad, b.grad)


def test_off_coarse_channel_bit_identical():
    seq, lg = rand_case(900, 2)
    n = 500
    rng = np.random.default_rng(3)
    seg = np.full(n, -100, dtype=np.int64)
    region = np.full(n, -100, dtype=np.int64)
    seg[:60] = 0                                         # exact background rows
    region[100:300] = rng.integers(0, 3, 200)           # coarse rows
    tg = dict(seg=torch.from_numpy(seg), region=torch.from_numpy(region))
    a = lg.clone().requires_grad_(True)
    r0 = coarse_channel_nll(seq, a, 128, 200, 700, tg, .1, cuts=(200, 700), checkpoint_bp=256)
    r0['loss'].backward()
    b = lg.clone().requires_grad_(True)
    r1 = GR.coarse_channel_nll(seq, b, 128, 200, 700, tg, .1, cuts=(200, 700), checkpoint_bp=256)
    r1['loss'].backward()
    assert r0['loss'].item() == r1['loss'].item() and torch.equal(a.grad, b.grad)


def test_off_decode_equals_frozen_skyline():
    from m29_skyline_reuse import decode_full_grammar_skyline
    for seed in (4, 5):
        seq, lg = rand_case(6000, seed, 3.)
        lg[:, 0] += 2.
        r0 = decode_full_grammar_skyline(seq, lg, 128, checkpoint_bp=1000)
        r1 = GR.decode_skyline(seq, lg, 128, checkpoint_bp=1000)
        assert r0['chains'] == r1['chains'] and r0['gain'] == r1['gain'] and r0['chain_gains'] == r1['chain_gains']


def test_disabled_states_refuse_bonus():
    seq, lg = rand_case(50, 6)
    with pytest.raises(ValueError):
        GR.log_partition(seq, lg, 8, P.bonus(TH, 1, P0), GR.OFF)


# ---- explicit enumeration ------------------------------------------------------------------------------------------------------------
def enumerate_genes(seq, logits, cap, intron_term):
    """Every complete gene (start, end, exons, score) of the frozen grammar, intron histories kept apart; score = exon gains (frozen
    `_TensorMargins`) + intron role gains + intron_term(k) per intron (k = length - 5)."""
    seq = seq.upper()
    L = len(seq)
    sc, orf = _TensorMargins(logits.detach()), ORFIndex(seq)
    starts = [i for i in range(L - 2) if seq[i:i + 3] == 'ATG']
    acceptors = [i + 2 for i in range(L - 2) if seq[i:i + 2] == 'AG']
    donors = {i for i in range(1, L - 1) if seq[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(L - 2) if seq[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    genes = []

    def extend(s, a, carry, count, score, exons):
        first = count == 0
        for b, terminal in events:
            if b <= a or ((first or terminal) and b - a < 3):
                continue
            rem = orf.extend(a, b, carry, terminal)
            if rem is None:
                continue
            v = score + float(sc.exon_gain(a, b, len(carry), first, terminal))
            ex = exons + [(a, b)]
            if terminal:
                genes.append((s, b, tuple(ex), v))
            elif count + 1 < cap:
                role = 4 + (len(rem) - 1) % 3
                for a2 in acceptors:
                    if a2 >= b + gap:
                        extend(s, a2, rem, count + 1, v + float(sc.prefix[role, a2] - sc.prefix[role, b]) + intron_term(a2 - b - gap), ex)
    for s in starts:
        extend(s, s, '', 0, 0., [])
    return genes


def set_value(genes, L, combine):
    """Sum (logsumexp) or max over all sets of nonoverlapping (adjacent allowed) genes, the empty set included."""
    by_start = {}
    for g in genes:
        by_start.setdefault(g[0], []).append(g)
    V = [0.] * (L + 1)
    for p in range(L - 1, -1, -1):
        vals = [V[p + 1]] + [g[3] + V[g[1]] for g in by_start.get(p, [])]
        V[p] = combine(vals)
    return V[0]


def lse(v):
    m = max(v)
    return m + math.log(sum(math.exp(x - m) for x in v))


def test_enumeration_matches_frozen_when_off():
    for seed in range(4):
        seq, lg = planted(seed)
        genes = enumerate_genes(seq, lg, 8, lambda k: 0.)
        z = set_value(genes, len(seq), lse)
        r = compiled_boundary_full_grammar_log_partition(seq, lg, 8, checkpoint_bp=16)
        assert len(genes) > 5
        assert abs(z - r['log_partition'].item()) < 1e-9


@pytest.mark.parametrize('il', [ON1, ON3])
def test_logz_matches_enumeration(il):
    r = il.short_stages
    bon = P.bonus(TH, r, P0)
    for seed in range(6):
        seq, lg = planted(seed)
        genes = enumerate_genes(seq, lg, 8, lambda k: P.log_len_sum(k, TH, r, P0))
        multi = sum(len(g[2]) > 2 for g in genes)
        z = set_value(genes, len(seq), lse)
        got = GR.log_partition(seq, lg, 8, bon, il, checkpoint_bp=7)['log_partition'].item()
        assert multi > 0 and abs(z - got) < 1e-9, (seed, z, got)


@pytest.mark.parametrize('il', [ON1, ON3])
def test_decode_matches_enumerated_max_and_is_legal(il):
    r = il.short_stages
    bon = P.bonus(TH, r, P0)
    for seed in range(6):
        seq, lg = planted(seed)
        lg = lg.float().double()
        lg[:, 0] -= 1.                                   # favour genes so the optimum is not the empty set
        genes = enumerate_genes(seq, lg, 8, lambda k: P.log_len_max(k, TH, r, P0))
        best = set_value(genes, len(seq), max)
        res = GR.decode_skyline(seq, lg, 8, bon, il, checkpoint_bp=9)
        assert abs(res['gain'] - best) < 1e-9, (seed, res['gain'], best)
        legal = {g[2]: g[3] for g in genes}
        for chain, gain in zip(res['chains'], res['chain_gains']):
            assert chain in legal and abs(legal[chain] - gain) < 1e-9
            for (a, b), (c, d) in zip(chain, chain[1:]):
                assert seq[b:b + 2] in ('GT', 'GC') and seq[c - 2:c] == 'AG' and c - b >= 5
        assert len(res['intron_branches']) == sum(len(c) - 1 for c in res['chains'])


def test_equal_rate_branches_reproduce_frozen_logz():
    """p_S = p_L = p0 (r = 1): the latent mixture sums to exactly the frozen grammar, for any pi."""
    lp0 = math.log(P0 / (1 - P0))
    seq, lg = rand_case(800, 7)
    z0 = GR.log_partition(seq, lg, 128)['log_partition'].item()
    for a in (-1., .3, 2.):
        th = torch.tensor([a, lp0, lp0], dtype=torch.float64)
        z1 = GR.log_partition(seq, lg, 128, P.bonus(th, 1, P0), ON1)['log_partition'].item()
        assert abs(z1 - z0) < 1e-9 * max(1, abs(z0))


# ---- gradients -------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('il', [ON1, ON3])
def test_gradcheck_logits_and_prior(il):
    seq, lg = planted(1)
    r = il.short_stages
    x = lg.double().clone().requires_grad_(True)
    th = TH.clone().requires_grad_(True)
    f = lambda x_, th_: GR.log_partition(seq, x_, 8, P.bonus(th_, r, P0), il, checkpoint_bp=5)['log_partition']
    assert torch.autograd.gradcheck(f, (x, th), eps=1e-6, atol=1e-6, rtol=1e-4)


def test_gradcheck_coarse_channel_prior():
    seq, lg = planted(2, 4)
    n = len(seq) - 10
    seg = np.full(n, -100, dtype=np.int64)
    region = np.full(n, -100, dtype=np.int64)
    region[5:25] = 1
    tg = dict(seg=torch.from_numpy(seg), region=torch.from_numpy(region))
    th = TH.clone().requires_grad_(True)
    x = lg.double().clone().requires_grad_(True)
    f = lambda x_, th_: GR.coarse_channel_nll(seq, x_, 8, 5, 5 + n, tg, .1, P.bonus(th_, 1, P0), ON1, checkpoint_bp=6)['loss']
    assert torch.autograd.gradcheck(f, (x, th), eps=1e-6, atol=1e-6, rtol=1e-4)
