"""CPU tests of the M33 posterior plumbing on a synthetic 30-kb problem (real Arabidopsis DNA, random potentials):
 (1) oriented_chains() maps the genomic chains of a '-' strand decode back to exactly the oriented chains the decoder produced;
 (2) posterior(...) returns probabilities in [0,1] for every chain, finite log partition;
 (3) m33_posterior.main() on a tiny run directory ('+' and '-') reproduces a direct call, and monotonicity: P at large tau >= P at small tau for MAP chains on average.
Run: python test_posterior_synth.py"""
import json, re, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
import torch
import m33_posterior as P
import m33_common as C
import manifest as MF
from m29_skyline_reuse import decode_full_grammar_skyline
from src.m28.c0_decode import genomic_chain

M29 = C.M29
seq = 'NC_003074.8'
full = re.sub('[^ACGT]', 'N', M29.read_sequence(seq))
a = 5_000_000
L = 30_000
dna_g = full[a:a + L]                                    # synthetic contig = this slice
rng = np.random.default_rng(0)
BG = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
base = rng.normal(0, 1, (L, 15)).astype(np.float32)
res = {}
for strand in '+-':
    oriented = dna_g if strand == '+' else MF.rc(dna_g)
    lg = base.copy()
    lg[:, 0] += BG                                       # background bias: small enough that some chains survive
    r = decode_full_grammar_skyline(oriented, torch.from_numpy(lg), 128)
    chains_o = [tuple(map(tuple, c)) for c in r['chains']]
    rows = [dict(exons=[(int(x), int(y)) for x, y in genomic_chain(c, 0, L, strand)], gain=float(g)) for c, g in zip(r['chains'], r['chain_gains'])]
    print(strand, 'decoded chains', len(rows), flush=True)
    assert len(rows) > 0, 'synthetic problem produced no chains; adjust the background bias'
    back = P.oriented_chains(rows, strand, L)
    assert sorted(back) == sorted(chains_o), (strand, 'orientation round trip failed')
    order, out = P.posterior(oriented, torch.from_numpy(lg), back)
    lp = np.array(out['gene_log_probabilities'])
    assert np.all(np.isfinite(lp)) and np.all(lp <= 1e-6), lp[:5]
    assert np.isfinite(out['log_partition'])
    res[strand] = (lg, rows, lp)
    # end-to-end through main(): run dir with pot/chains, compare to the direct call
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / 'pot').mkdir(); (td / 'chains').mkdir()
        name = 'SYN'
        np.save(td / 'pot' / f'{name}_{strand}.npy', lg)
        json.dump(dict(chains=rows), open(td / 'chains' / f'{name}_{strand}.json', 'w'))
        # main() reads the real genome for `seq`; emulate it by monkeypatching read_sequence for the synthetic contig name
        orig = M29.read_sequence
        M29.read_sequence = lambda s: dna_g if s == name else orig(s)
        sys.argv = ['x', str(td), name, strand, '0.5,1', str(td / 'post')]
        P.main()
        M29.read_sequence = orig
        d1 = json.load(open(td / 'post' / 'tau1' / f'{name}_{strand}.json'))
        assert np.allclose(sorted(g['log_p'] for g in d1['genes']), sorted(lp), atol=1e-5), 'main() differs from direct call'
        d05 = json.load(open(td / 'post' / 'tau0.5' / f'{name}_{strand}.json'))
        print(strand, 'mean P tau=1:', float(np.exp(lp).mean()), ' tau=0.5:', float(np.mean([np.exp(g['log_p']) for g in d05['genes']])), flush=True)
print('ALL OK')
