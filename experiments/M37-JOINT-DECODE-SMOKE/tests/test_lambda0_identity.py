"""lambda = 0 joint decoding == A-alone decoding, chain by chain.
 (a) synthetic: random sequence + random potentials, decoded directly and after blend(lambda=0) with a random projection -> identical chains
     (frozen decoder `decode_full_grammar_skyline`, exactly the call used in m33_infer / m37_joint);
 (b) real: outputs/M37-JOINT-DECODE-SMOKE/joint/lambda0/chains/NC_089044.1_{+,-}.json vs M33 val chains (skipped when absent)."""
import json, random, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import m37_common as C


def test_synthetic():
    import torch
    from m29_skyline_reuse import decode_full_grammar_skyline
    from m37_joint import blend
    rng = random.Random(5); nr = np.random.default_rng(5)
    s = ''.join(rng.choice('ACGT') for _ in range(20_000))
    pot = nr.standard_normal((len(s), 15)).astype(np.float32) * 2
    proj = nr.standard_normal((len(s), 15)).astype(np.float32)
    mask = nr.random(len(s)) < .7
    z = blend(pot, proj, mask, 0.)
    assert np.array_equal(z, pot)
    r1 = decode_full_grammar_skyline(s, torch.from_numpy(pot), 128)
    r2 = decode_full_grammar_skyline(s, torch.from_numpy(z), 128)
    assert r1['chains'] == r2['chains'] and r1['gain'] == r2['gain'], (len(r1['chains']), len(r2['chains']))
    print('test_synthetic OK', len(r1['chains']), 'chains')


def test_real():
    d = C.OUT / 'joint' / 'lambda0' / 'chains'
    if not all((d / f'{C.A_CONTIG}_{st}.json').exists() for st in '+-'):
        print('SKIP test_real: lambda0 chains absent'); return
    for st in '+-':
        j, _ = C.load_chain_json(d / f'{C.A_CONTIG}_{st}.json')
        a, _ = C.load_chain_json(C.A_CHAINS[st])
        assert [x[0] for x in j] == [x[0] for x in a], (st, len(j), len(a))
        assert np.allclose([x[1] for x in j], [x[1] for x in a])
        print('test_real OK', st, len(j), 'chains identical to M33 val')


if __name__ == '__main__':
    test_synthetic(); test_real()
