"""Parallel window DP (partition-function loss + gradient wrt the 15-role potentials), M29-A1.

Per window:  loss = log Z_all - log Z_obs  (frozen `coarse_channel_nll`, epsilon .1, 128 exons, 256-bp checkpoints) and
d loss / d logits (FP32, from the frozen compiled reverse pass).  Only change vs. the frozen path: the M29 reusable alpha arena in
`_reverse_messages` (bit-identical gradients, tested) and process-level parallelism across the windows of a batch.
"""
import hashlib
import multiprocessing as mp
import os
import time
import numpy as np


def _install():
    import torch
    torch.set_num_threads(1)
    import a1_common  # noqa: F401
    import src.m28.full_grammar_boundary_compiled as BC
    from m29_reverse_arena import _reverse_messages_arena
    BC._reverse_messages = _reverse_messages_arena
    return torch, BC


def window_loss_grad(context, logits, core_start, core_end, targets, epsilon=.1, arena=True):
    """context: str; logits: float32 numpy [L,15] (head + fixed offset); targets: dict(seg,region) numpy.
    Returns dict(loss, grad float32 [L,15], seconds)."""
    import torch
    from src.m28.full_grammar_coarse_channel import coarse_channel_nll
    t0 = time.perf_counter()
    L, n = len(context), core_end - core_start
    cuts = tuple(x for x in (core_start, core_end) if 0 < x < L)
    lg = torch.from_numpy(np.ascontiguousarray(logits, dtype=np.float32)).clone().requires_grad_(True)
    tg = {k: torch.as_tensor(v) for k, v in targets.items()}
    res = coarse_channel_nll(context, lg, 128, core_start, core_end, tg, epsilon, cuts=cuts, checkpoint_bp=256)
    res['loss'].backward()
    return dict(loss=float(res['loss'].detach()), grad=lg.grad.numpy().copy(), seconds=time.perf_counter() - t0,
                log_Z_all=float(res['log_partition_all'].detach()), log_Z_obs=float(res['log_partition_observed'].detach()))


_CANARY = None


def _canary():
    """Tiny deterministic DP (1.2 kb); a digest change in a live worker means that core corrupts FP64/interpreter state."""
    import torch
    from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
    rng = np.random.default_rng(99)
    seq = ''.join(rng.choice(list('ACGT'), 1200))
    lg = torch.from_numpy((rng.standard_normal((1200, 15)) * 2).astype(np.float32)).requires_grad_(True)
    r = compiled_boundary_full_grammar_log_partition(seq, lg, 128, checkpoint_bp=256)
    r['log_partition'].backward()
    return hashlib.sha256(np.float64(r['log_partition'].item()).tobytes() + lg.grad.numpy().tobytes()).hexdigest()


def _init(cpus):
    global _CANARY
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    _install()
    _CANARY = _canary()


def _job(args):
    context, logits, cs, ce, targets, epsilon = args
    r = window_loss_grad(context, logits, cs, ce, targets, epsilon)
    c = _canary()
    if c != _CANARY:
        raise RuntimeError(f'worker canary digest changed (pid {os.getpid()}, cpu {os.sched_getaffinity(0)}): possible core fault')
    r['pid'] = os.getpid()
    return r


class DPPool:
    def __init__(self, workers, cpus=None):
        self.pool = mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(cpus,))

    def submit(self, context, logits, core_start, core_end, targets, epsilon=.1):
        return self.pool.apply_async(_job, ((context, logits, core_start, core_end, targets, epsilon),))

    def close(self):
        self.pool.terminate()
        self.pool.join()
