"""ont stability self-check (run before any long task).
Deterministic CPU numpy/numba workloads in N parallel processes, repeated; every repetition must reproduce the same
SHA-256.  Also a deterministic GPU matmul check and a frozen-DP (numba FP64) repeatability check on a fixed random input.
Usage: a1_selfcheck.py <minutes> <procs>
"""
import hashlib, json, multiprocessing as mp, os, sys, time
import numpy as np


def cpu_job(args):
    seed, seconds = args
    from numba import njit

    @njit(cache=False)
    def kernel(a, n):
        s = 0.0
        for i in range(n):
            s += np.sin(a[i % len(a)]) * np.exp(-a[(i * 7) % len(a)] * 1e-3)
        return s

    rng = np.random.default_rng(seed)
    a = rng.standard_normal(4096)
    m = rng.standard_normal((512, 512))
    ref = None
    reps, bad = 0, 0
    t0 = time.time()
    while time.time() - t0 < seconds or reps < 3:
        h = hashlib.sha256()
        h.update(np.float64(kernel(a, 2_000_000)).tobytes())
        x = m.copy()
        for _ in range(4):
            x = x @ m
            x /= np.abs(x).max()
        h.update(x.tobytes())
        h.update(np.sort(rng.standard_normal(1)).tobytes() if False else b'')
        d = h.hexdigest()
        if ref is None:
            ref = d
        elif d != ref:
            bad += 1
        reps += 1
    return dict(seed=seed, reps=reps, mismatches=bad, digest=ref)


def gpu_check():
    import torch
    if not torch.cuda.is_available():
        return dict(gpu=False)
    g = torch.Generator(device='cuda').manual_seed(0)
    a = torch.randn(4096, 4096, device='cuda', generator=g)
    ref = None
    bad = 0
    for _ in range(200):
        c = (a @ a).double().sum().item()
        if ref is None:
            ref = c
        elif c != ref:
            bad += 1
    return dict(gpu=True, reps=200, mismatches=bad, value=ref)


def dp_check():
    sys.path[:0] = [os.path.expanduser('~/coding-rna-reorg/pipeline/m28'), os.path.expanduser('~/coding-rna-reorg/pipeline/m28/scripts')]
    import torch
    from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
    rng = np.random.default_rng(1)
    seq = ''.join(rng.choice(list('ACGT'), 6000))
    lg = torch.from_numpy((rng.standard_normal((6000, 15)) * 2).astype(np.float32)).requires_grad_(True)
    ref, bad = None, 0
    for _ in range(6):
        lg.grad = None
        r = compiled_boundary_full_grammar_log_partition(seq, lg, 128, checkpoint_bp=256)
        r['log_partition'].backward()
        d = hashlib.sha256(np.float64(r['log_partition'].item()).tobytes() + lg.grad.numpy().tobytes()).hexdigest()
        if ref is None:
            ref = d
        elif d != ref:
            bad += 1
    return dict(reps=6, mismatches=bad, digest=ref)


if __name__ == '__main__':
    minutes, procs = float(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    # same seed for pairs of processes -> digests must also agree across processes
    with mp.get_context('spawn').Pool(procs) as pool:
        res = pool.map(cpu_job, [(i % max(1, procs // 2), minutes * 60) for i in range(procs)])
    by_seed = {}
    for r in res:
        by_seed.setdefault(r['seed'], set()).add(r['digest'])
    out = dict(cpu=res, cross_process_digest_disagreements=sum(len(v) > 1 for v in by_seed.values()),
               gpu=gpu_check(), dp=dp_check(), seconds=time.time() - t0)
    print(json.dumps(out, indent=1))
