"""CPU unit test of the one-hot floor backend (protocol 32 arm L, `M39_BACKEND=onehot`; user 2026-11-12).  No GPU, no language model.
Checks
  1 features : `onehot_tokens` equals a naive per-base construction (A,C,G,T -> dims 4j..4j+3 of token t for oriented base 6t+j; N / lowercase-free
               padding -> 0; dims >= 24 are 0) on random strings with N, full and short windows.
  2 backend  : load(lora=False) gives no encoder and no parameters; load(lora=True) is refused; encoder_cfg is JSON-able; init = M32 runA head.
  3 draws    : real draws of schedule_I10_k50U (both strands) -> pc2 geometry -> window_pots on the CPU with the M32 runA head: finite [len, 15] per
               window; the loss gradient reaches the head projection columns 0..23 and only those input columns.
  4 regions  : region_potentials on a 20-kb slice of a new-species validation crop (both strands): finite [L, 15].
Output: JSON on stdout (and --out FILE).
usage: M39_BACKEND=onehot ~/envs/plantcad2/bin/python m39_oh_cpu_check.py [--out FILE]
"""
import argparse
import json
import os
import random

os.environ.setdefault('M39_BACKEND', 'onehot')
assert os.environ['M39_BACKEND'] == 'onehot', 'run with M39_BACKEND=onehot'
import torch                                                  # noqa: E402

import m39_train_common as T                                  # noqa: E402
import m39_backend as BK                                      # noqa: E402

A, G, MF, C, MD = T.A, T.G, T.MF, T.C, T.MD


def naive(s):
    f = torch.zeros(A.HEAD_TOKENS, A.HIDDEN)
    for p, c in enumerate(s):
        if c in 'ACGT':
            f[p // 6, 4 * (p % 6) + 'ACGT'.index(c)] = 1
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out')
    args = ap.parse_args()
    rng = random.Random(0)
    out = dict(backend=T.BACKEND, init=str(T.INIT), problems=[])
    # 1 features
    strings = [''.join(rng.choice('ACGTN' if k % 3 else 'ACGT') for _ in range(n)) for k, n in enumerate((8192, 8192, 6001, 1003, 7, 1))]
    strings.append('N' * 500)
    bad = 0
    for s, f in zip(strings, BK.onehot_tokens(strings, 'cpu')):
        if f.dtype != torch.bfloat16 or tuple(f.shape) != (A.HEAD_TOKENS, A.HIDDEN) or not torch.equal(f.float(), naive(s)):
            bad += 1
    out['features'] = dict(strings=len(strings), mismatches=bad, dims=BK.OH_DIMS)
    if bad:
        out['problems'].append('onehot_tokens mismatch')
    # 2 backend
    init = torch.load(T.INIT, map_location='cpu', weights_only=False)
    BK.check_init(init)
    cfg = dict(seed=0, **BK.encoder_cfg(init['config']))
    json.dumps(cfg)
    ctx, params = BK.load(cfg, lora=False, device='cpu')
    try:
        BK.load(cfg, lora=True, device='cpu')
        out['problems'].append('load(lora=True) accepted')
    except AssertionError:
        pass
    if ctx.enc is not None or params or BK.lora_state(ctx):
        out['problems'].append('onehot ctx carries an encoder / parameters')
    out['encoder_cfg'] = cfg
    head = BK.new_head(cfg, 0)
    head.load_state_dict(init['model'])
    head.train()
    # 3 draws
    rows = T.read_schedule(T.SCHEDULES / 'schedule_I10_k50U.jsonl')
    picks = [r for r in rows[:400] if r['strand'] == '+'][:2] + [r for r in rows[:400] if r['strand'] == '-'][:2]
    tab = C.species_table()
    proj = [p for n, p in head.named_parameters() if p.dim() == 2 and p.shape[1] == A.HIDDEN]
    assert proj, 'no head parameter with input width HIDDEN'
    draws = []
    for r in picks:
        sp = tab[r['species']]
        chrom = MF.read_allowed(sp['dir'] / 'genome.fa', {r['seqid']: sp['lengths'][r['seqid']]})[r['seqid']]
        geo = G.context_geometry(r['start'], r['end'], len(chrom))
        pg = BK.geometry(*geo['genomic_context'])
        it = dict(strand=r['strand'], **BK.item_windows(chrom, pg, r['strand']))
        assert set(it) == {'strand', 'seqs'}
        head.zero_grad()
        pots = BK.window_pots(head, ctx, it, list(range(len(it['seqs']))), True, device='cpu')
        ok = all(tuple(x.shape) == (len(s), 15) and bool(torch.isfinite(x).all()) for x, s in zip(pots, it['seqs']))
        sum(x.float().pow(2).mean() for x in pots).backward()
        g = proj[0].grad.abs().sum(0)
        grad_ok = bool((g[:BK.OH_DIMS] > 0).any()) and float(g[BK.OH_DIMS:].sum()) == 0.0
        draws.append(dict(id=r['id'], species=r['species'], strand=r['strand'], windows=len(pots), finite_shapes=ok, grad_only_first_24=grad_ok))
        if not (ok and grad_ok):
            out['problems'].append(dict(draw=r['id'], finite_shapes=ok, grad=grad_ok))
    out['draws'] = draws
    # 4 regions
    head.eval()
    MC = T.register()
    contig = T.val_contig(T.new_keys()[0])
    sl = MC.read_sequence(contig)[:20000]
    import re
    sl = re.sub('[^ACGT]', 'N', sl)
    pots, n = BK.region_potentials(head, ctx, sl, '+-', batch=4, device='cpu')
    reg_ok = all(tuple(pots[s].shape) == (len(sl), 15) and bool(torch.isfinite(pots[s]).all()) for s in '+-')
    out['region'] = dict(contig=contig, length=len(sl), windows=n, finite_shapes=reg_ok)
    if not reg_ok:
        out['problems'].append('region_potentials')
    out['ok'] = not out['problems']
    s = json.dumps(out, indent=1, default=str)
    if args.out:
        open(args.out, 'w').write(s)
    print(s)


if __name__ == '__main__':
    main()
