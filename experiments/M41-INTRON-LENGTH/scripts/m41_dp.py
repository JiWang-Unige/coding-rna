"""M41 parallel window DP = M29 `a1_dp` (partition-function loss + gradient wrt the 15-role potentials, epsilon .1, 128 exons, 256-bp checkpoints, process
pool, canary) on the M41 automaton, also returning the gradient wrt the intron-length bonus vector.  `a1_dp` is imported, not edited; its frozen canary
(the frozen automaton) still guards every worker core."""
import os
import time

import numpy as np

import a1_dp as DP0                                  # M29-A1 (unchanged)


def window_loss_grad(context, logits, core_start, core_end, targets, epsilon, bonus, short_stages):
    """context: str; logits float32 [L,15]; targets dict(seg, region); bonus float64 [N_TYPES]; short_stages 0 = states off.
    -> dict(loss, grad float32 [L,15], grad_bonus float64 [N_TYPES], seconds, log_Z_all, log_Z_obs)"""
    import torch
    import m41_grammar as GR
    t0 = time.perf_counter()
    L = len(context)
    il = GR.IntronLength(True, short_stages) if short_stages else GR.OFF
    cuts = tuple(x for x in (core_start, core_end) if 0 < x < L)
    lg = torch.from_numpy(np.ascontiguousarray(logits, dtype=np.float32)).clone().requires_grad_(True)
    bo = torch.from_numpy(np.asarray(bonus, dtype=np.float64)).clone().requires_grad_(True)
    tg = {k: torch.as_tensor(v) for k, v in targets.items()}
    res = GR.coarse_channel_nll(context, lg, 128, core_start, core_end, tg, epsilon, bo, il, cuts=cuts, checkpoint_bp=256)
    res['loss'].backward()
    return dict(loss=float(res['loss'].detach()), grad=lg.grad.numpy().copy(), grad_bonus=bo.grad.numpy().copy(), seconds=time.perf_counter() - t0,
                log_Z_all=float(res['log_partition_all'].detach()), log_Z_obs=float(res['log_partition_observed'].detach()))


def _init(cpus):
    DP0._init(cpus)                                  # affinity, frozen-path install, canary digest


def _job(args):
    context, logits, cs, ce, targets, epsilon, bonus, r = args
    out = window_loss_grad(context, logits, cs, ce, targets, epsilon, bonus, r)
    c = DP0._canary()
    if c != DP0._CANARY:
        raise RuntimeError(f'worker canary digest changed (pid {os.getpid()}, cpu {os.sched_getaffinity(0)}): possible core fault')
    out['pid'] = os.getpid()
    return out


class DPPool:
    def __init__(self, workers, cpus=None):
        import multiprocessing as mp
        self.pool = mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(cpus,))

    def submit(self, context, logits, core_start, core_end, targets, epsilon, bonus, short_stages):
        return self.pool.apply_async(_job, ((context, logits, core_start, core_end, targets, epsilon, bonus, short_stages),))

    def close(self):
        self.pool.terminate()
        self.pool.join()
