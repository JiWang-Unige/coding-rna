"""Parallel (process-pool, arena-patched) window DP == frozen serial implementation on real TRAIN windows with PlantCAD2-sized contexts
(random 15-role logits, no offset): loss and gradient must agree bitwise.  Same method as the M29-A1 test."""
import numpy as np
import torch
import m32_common as M
import a1_dp as DP
import src.m28.full_grammar_boundary_compiled as BC
from src.m28.full_grammar_coarse_channel import coarse_channel_nll


def serial_frozen(context, logits, cs, ce, targets):
    assert BC._reverse_messages.__name__ == '_reverse_messages'
    torch.set_num_threads(1)
    lg = torch.from_numpy(logits).clone().requires_grad_(True)
    cuts = tuple(x for x in (cs, ce) if 0 < x < len(context))
    tg = {k: torch.as_tensor(v) for k, v in targets.items()}
    r = coarse_channel_nll(context, lg, 128, cs, ce, tg, .1, cuts=cuts, checkpoint_bp=256)
    r['loss'].backward()
    return float(r['loss'].detach()), lg.grad.numpy().copy()


def test_parallel_equals_serial_and_heldout_never_sampled():
    windows, strata, catalogs, original = M.PP.load_population()
    sequences, _ = M.C.P.read_train_DNA()
    sched = M.schedule(4608 * 8, strata, original)
    assert not any(windows[d['window_id']]['seqid'] in M.HELDOUT for d in sched)
    rng = np.random.default_rng(5)
    jobs, quota = [], {}
    for d in sched[:400]:
        w = windows[d['window_id']]
        key = (w['species'], w['strand'])
        if quota.get(key, 0) >= 1:
            continue
        chromosome, geometry, context, offset, targets, counts = M.example(w, sequences, catalogs)
        if counts[0] < 300:
            continue
        quota[key] = 1
        head = (rng.standard_normal((len(context), 15)) * 1.5).astype(np.float32)
        jobs.append((context, head, offset, offset + w['valid_bases'], {k: np.asarray(targets[k]) for k in ('seg', 'region')}))
    assert len(jobs) == 4
    pool = DP.DPPool(4)
    par = [f.get() for f in [pool.submit(*j) for j in jobs]]
    pool.close()
    for j, p in zip(jobs, par):
        loss, grad = serial_frozen(*j)
        assert p['loss'] == loss
        assert np.array_equal(p['grad'], grad)
        assert p['loss'] > 0 and np.abs(p['grad']).max() > 0
    print('losses', [round(p['loss'], 2) for p in par], 'dp_s', [round(p['seconds'], 1) for p in par])
