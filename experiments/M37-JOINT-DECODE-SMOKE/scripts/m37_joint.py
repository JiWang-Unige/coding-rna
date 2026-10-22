"""Potential projection B -> A and potential-averaged joint decoding of A (protocol §4.2).
pot_A'(lambda) = (1 - lambda) * pot_A + lambda * proj(pot_B) on mapped A bases; pot_A elsewhere.  Oriented 15-role potentials: the '+' file of a
contig is indexed by forward position, the '-' file by reverse-complement position (L-1-i).  Projection of an A base a mapped to B base b:
  same strand   : proj+[a] = potB+[b]            proj-[LA-1-a] = potB-[LB-1-b]
  inverted block: proj+[a] = potB-[LB-1-b]       proj-[LA-1-a] = potB+[b]        (strand potentials swapped, positions reversed)
Decoding: frozen exact whole-contig skyline (`decode_full_grammar_skyline`, 128 exons), exactly the decode call of m33_infer / M31 c1_infer.
lambda = 0 must reproduce the A-alone chains of M33 (tests/test_lambda0_identity.py).
Usage: m37_joint.py --map OUT/align/map.npz --out OUT/joint --lambdas 0,0.25,0.5 [--workers 4 --cpus 25-31 --keep-pot]
"""
import argparse, json, multiprocessing as mp, os, re, sys, time
from pathlib import Path
import numpy as np
import torch
import m37_common as C
import m31_common as M31
import manifest as MF
C.register()


def project(potB_p, potB_m, a2b, a2b_s):
    LA, LB = len(a2b), len(potB_p)
    mapped = a2b >= 0
    proj_p = np.full((LA, 15), np.nan, np.float32); proj_m = np.full((LA, 15), np.nan, np.float32)
    ia = np.nonzero(mapped & (a2b_s > 0))[0]; ib = a2b[ia]
    proj_p[ia] = potB_p[ib]; proj_m[LA - 1 - ia] = potB_m[LB - 1 - ib]
    ia = np.nonzero(mapped & (a2b_s < 0))[0]; ib = a2b[ia]
    proj_p[ia] = potB_m[LB - 1 - ib]; proj_m[LA - 1 - ia] = potB_p[ib]
    return {'+': proj_p, '-': proj_m}, {'+': mapped, '-': mapped[::-1].copy()}


def blend(potA, proj, mask, lam):
    out = potA.copy()
    out[mask] = (1. - lam) * potA[mask] + lam * proj[mask]
    assert np.isfinite(out).all()
    return out


def _init(cpus):
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(1)
    import a1_common
    a1_common.repoint()
    C.register()


def decode_job(seq, strand, pot_path, out_json, keep):
    from m29_skyline_reuse import decode_full_grammar_skyline
    from src.m28.c0_decode import genomic_chain
    M = C.register()
    lg = np.load(pot_path)
    full = re.sub('[^ACGT]', 'N', M.read_sequence(seq))
    L = len(full)
    assert lg.shape == (L, 15) and lg.dtype == np.float32
    oriented = full if strand == '+' else MF.rc(full)
    t0 = time.perf_counter()
    res = decode_full_grammar_skyline(oriented, torch.from_numpy(lg), 128)
    sec = time.perf_counter() - t0
    chains = [dict(exons=[(int(x), int(y)) for x, y in genomic_chain(ch, 0, L, strand)], gain=float(g)) for ch, g in zip(res['chains'], res['chain_gains'])]
    json.dump(dict(seq=seq, start=0, end=L, strand=strand, chains=chains, total_gain=res['gain'], decode_seconds=sec), open(out_json, 'w'))
    if not keep:
        os.remove(pot_path)
    return out_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--map', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--lambdas', default='0,0.25,0.5'); ap.add_argument('--workers', type=int, default=4); ap.add_argument('--cpus', default='')
    ap.add_argument('--keep-pot', action='store_true', help='keep the blended potentials (default: delete after decoding; recomputable from A/B pots + map)')
    ap.add_argument('--potB', default=str(C.OUT / 'inferB' / 'pot'))
    a = ap.parse_args()
    cpus = M31.parse_cpus(a.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    lams = [float(x) for x in a.lambdas.split(',')]
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    m = np.load(a.map)
    a2b, a2b_s = m['a2b'], m['a2b_strand']
    potA = {st: np.load(C.A_POT[st]) for st in '+-'}
    potB = {st: np.load(Path(a.potB) / f'{C.B_CONTIG}_{st}.npy') for st in '+-'}
    assert potA['+'].shape[0] == len(a2b) and potB['+'].shape[0] == len(m['b2a'])
    proj, mask = project(potB['+'], potB['-'], a2b, a2b_s)
    del potB
    pool = mp.get_context('spawn').Pool(a.workers, initializer=_init, initargs=(cpus,))
    jobs, log = [], []
    for lam in lams:
        d = out / f'lambda{lam:g}'
        (d / 'pot').mkdir(parents=True, exist_ok=True); (d / 'chains').mkdir(exist_ok=True)
        for st in '+-':
            cj = d / 'chains' / f'{C.A_CONTIG}_{st}.json'
            if cj.exists():
                continue
            pp = d / 'pot' / f'{C.A_CONTIG}_{st}.npy'
            if not pp.exists():
                b = blend(potA[st], proj[st], mask[st], lam)
                np.save(str(pp) + '.tmp.npy', b); os.replace(str(pp) + '.tmp.npy', pp); del b
            jobs.append((lam, st, pool.apply_async(decode_job, (C.A_CONTIG, st, str(pp), str(cj), a.keep_pot))))
            print(json.dumps(dict(event='submitted', **{'lambda': lam}, strand=st)), flush=True)
    del potA, proj
    for lam, st, j in jobs:
        f = j.get()
        log.append(dict(**{'lambda': lam}, strand=st, chains=f, decode_seconds=json.load(open(f))['decode_seconds']))
        print(json.dumps(log[-1]), flush=True)
    pool.close(); pool.join()
    json.dump(dict(lambdas=lams, map=a.map, mapped_A_bases=int(mask['+'].sum()), LA=len(a2b), decodes=log), open(out / 'DONE.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
