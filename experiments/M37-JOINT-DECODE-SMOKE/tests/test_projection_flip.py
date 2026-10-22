"""Projection of oriented potentials through the base map (m37_joint.project / blend):
 (1) formulas for same-strand and inverted blocks, unmapped bases untouched, lambda = 0 -> pot_A exactly, lambda = 1 -> projected values;
 (2) semantic check: potentials defined as a function of the oriented local sequence; with B = rc(A) (one whole inverted block) the projection
     must reproduce A's own oriented potentials on both strands, i.e. the strand swap + position reversal is the right convention."""
import random, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
COMP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(COMP)[::-1]


def seq_pot(s):
    """15-role 'potential' that is a deterministic function of the oriented 5-mer around each position"""
    L = len(s); out = np.zeros((L, 15), np.float32)
    code = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 0}
    for i in range(L):
        h = 0
        for j in range(i - 2, i + 3):
            h = h * 5 + (code[s[j]] if 0 <= j < L else 0)
        rng = np.random.default_rng(h)
        out[i] = rng.standard_normal(15)
    return out


def test_formulas():
    from m37_joint import project, blend
    rng = np.random.default_rng(1)
    LA, LB = 50, 60
    a2b = np.full(LA, -1); a2b_s = np.zeros(LA, np.int8)
    a2b[5:15] = np.arange(20, 30); a2b_s[5:15] = 1                       # same strand block
    a2b[30:40] = np.arange(55, 45, -1); a2b_s[30:40] = -1                # inverted block: a=30 <-> b=55, a=39 <-> b=46
    pB = {st: rng.standard_normal((LB, 15)).astype(np.float32) for st in '+-'}
    pA = {st: rng.standard_normal((LA, 15)).astype(np.float32) for st in '+-'}
    proj, mask = project(pB['+'], pB['-'], a2b, a2b_s)
    for a in range(5, 15):
        b = a2b[a]
        assert np.array_equal(proj['+'][a], pB['+'][b]) and np.array_equal(proj['-'][LA - 1 - a], pB['-'][LB - 1 - b])
    for a in range(30, 40):
        b = a2b[a]
        assert np.array_equal(proj['+'][a], pB['-'][LB - 1 - b]) and np.array_equal(proj['-'][LA - 1 - a], pB['+'][b])
    assert np.isnan(proj['+'][0]).all() and np.isnan(proj['-'][LA - 1]).all()
    assert mask['+'].sum() == 20 and np.array_equal(mask['-'], mask['+'][::-1])
    for st in '+-':
        z = blend(pA[st], proj[st], mask[st], 0.)
        assert np.array_equal(z, pA[st])                                  # lambda = 0: bit-identical
        o = blend(pA[st], proj[st], mask[st], 1.)
        assert np.array_equal(o[mask[st]], proj[st][mask[st]]) and np.array_equal(o[~mask[st]], pA[st][~mask[st]])
        h = blend(pA[st], proj[st], mask[st], .5)
        assert np.allclose(h[mask[st]], .5 * pA[st][mask[st]] + .5 * proj[st][mask[st]]) and np.array_equal(h[~mask[st]], pA[st][~mask[st]])
    print('test_formulas OK')


def test_semantic_inversion():
    from m37_joint import project
    rng = random.Random(3)
    A = ''.join(rng.choice('ACGT') for _ in range(300)); B = rc(A)
    LA = LB = len(A)
    a2b = LA - 1 - np.arange(LA); a2b_s = np.full(LA, -1, np.int8)
    pA = {'+': seq_pot(A), '-': seq_pot(rc(A))}
    pB = {'+': seq_pot(B), '-': seq_pot(rc(B))}
    proj, mask = project(pB['+'], pB['-'], a2b, a2b_s)
    assert mask['+'].all()
    for st in '+-':
        assert np.array_equal(proj[st], pA[st]), st
    # same-strand identity map as the control
    a2b = np.arange(LA); a2b_s = np.ones(LA, np.int8)
    proj2, _ = project(pA['+'], pA['-'], a2b, a2b_s)
    for st in '+-':
        assert np.array_equal(proj2[st], pA[st])
    print('test_semantic_inversion OK')


if __name__ == '__main__':
    test_formulas(); test_semantic_inversion()
