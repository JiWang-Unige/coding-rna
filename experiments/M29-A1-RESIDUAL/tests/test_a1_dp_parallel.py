"""Parallel (process-pool, arena-patched) window DP == frozen serial implementation on >= 8 real TRAIN windows:
loss and gradient must agree (<= 1e-10 relative; in practice bitwise)."""
import json
import numpy as np
import pytest
import torch
import a1_common as K
import a1_dp as DP
import a1_model as MD     # noqa: F401  (repoint)
import a1_offset as O
import prepare_full_grammar_fresh_draws as PP
import fit_full_grammar_objective_contrast as C
import a1_train as T
import src.m28.full_grammar_boundary_compiled as BC
from src.m28.full_grammar_coarse_channel import coarse_channel_nll

N_WIN = 8


def serial_frozen(context, logits, cs, ce, targets):
    assert BC._reverse_messages.__name__ == '_reverse_messages', 'serial reference must use the frozen reverse pass'
    torch.set_num_threads(1)
    lg = torch.from_numpy(logits).clone().requires_grad_(True)
    L = len(context)
    cuts = tuple(x for x in (cs, ce) if 0 < x < L)
    tg = {k: torch.as_tensor(v) for k, v in targets.items()}
    r = coarse_channel_nll(context, lg, 128, cs, ce, tg, .1, cuts=cuts, checkpoint_bp=256)
    r['loss'].backward()
    return float(r['loss'].detach()), lg.grad.numpy().copy()


def test_parallel_equals_serial():
    plan = [json.loads(x) for x in (K.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    tracks = {s: O.load_tracks(s) for s in {w['seqid'] for w in windows.values()}}
    # 8 windows whose core contains gene roles (2 per species/strand) + 1 fully-unknown window (loss and gradient must be exactly 0).
    # Gene-free windows also give loss 0.0 in FP64 under the Carbon-A offset (genes outside genes are suppressed by ~exp(-1e3)).
    quota, idx = {}, []
    zero = next((i for i, r in enumerate(plan) if r['fine_positions'] + r['coarse_positions'] == 0), None)
    for i, r in enumerate(plan[:600]):
        w = windows[r['window_id']]
        key = (w['species'], w['strand'])
        if quota.get(key, 0) < 2 and r['fine_positions'] > 200:
            chromosome, geometry, context, offset, targets, counts = T.example(w, sequences, catalogs)
            if int((torch.as_tensor(targets['seg'])[:w['valid_bases']] >= 1).sum()) > 300:
                quota[key] = quota.get(key, 0) + 1
                idx.append(i)
        if len(idx) == N_WIN and zero is not None:
            break
    assert len(idx) == 8, idx
    assert zero is not None, 'plan has no all-unknown window'
    idx = idx + [zero]
    jobs = []
    rng = np.random.default_rng(5)
    for i in idx:
        w = windows[plan[i]['window_id']]
        chromosome, geometry, context, offset, targets, counts = T.example(w, sequences, catalogs)
        a, b = geometry['genomic_context']
        off = O.context_offset(tracks[w['seqid']], w['strand'], a, b, len(chromosome))
        head = (rng.standard_normal(off.shape) * 3.).astype(np.float32)       # a large random residual: under the Carbon-A offset the label path is already ~1.0 (loss 0.0 in FP64) for most gene windows
        jobs.append((context, (off + head).astype(np.float32), offset, offset + w['valid_bases'],
                     {k: np.asarray(targets[k]) for k in ('seg', 'region')}))
    pool = DP.DPPool(9)
    futs = [pool.submit(*j) for j in jobs]
    par = [f.get() for f in futs]
    pool.close()
    assert par[-1]['loss'] == 0.0 and not np.any(par[-1]['grad'])
    for j, p in zip(jobs[:-1], par[:-1]):
        pass
        loss, grad = serial_frozen(*j)
        assert p['loss'] == loss or abs(p['loss'] - loss) <= 1e-10 * abs(loss)
        denom = np.abs(grad).max()
        assert np.abs(p['grad'] - grad).max() <= 1e-10 * max(denom, 1e-30)
        assert np.array_equal(p['grad'], grad)       # bitwise in practice
    assert sum(p['loss'] > 0 and np.abs(p['grad']).max() > 0 for p in par[:-1]) >= 6
    print('windows', len(jobs), 'losses', [round(p['loss'], 3) for p in par], 'dp_s', [round(p['seconds'], 1) for p in par])
