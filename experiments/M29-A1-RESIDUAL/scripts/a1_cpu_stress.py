"""Per-logical-CPU pinned determinism stress (diagnose a suspected faulty core on ont).
dmesg shows ~30 userspace segfaults 'likely on CPU 0/1 (core 0)' since uptime ~704,000 s; this pins one worker to each
logical CPU, runs interpreter-heavy + numpy + numba workloads with a fixed seed and checks the digest every repetition.
Usage: a1_cpu_stress.py <seconds> [cpu,cpu,...]   (default: all CPUs)
"""
import hashlib, json, multiprocessing as mp, os, sys, time, traceback
from collections import Counter


def py_workload(seed):
    # interpreter-heavy: dict/set/list/str/Counter/sorting/closures, fully deterministic
    x = seed
    d, s, l = {}, set(), []
    for i in range(60000):
        x = (x * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        k = x >> 40
        d[k] = d.get(k, 0) + 1
        s.add(k % 977)
        l.append(str(k)[::-1])
    c = Counter(l)
    top = sorted(c.items(), key=lambda t: (-t[1], t[0]))[:50]
    f = (lambda a: (lambda b: a + b))
    acc = 0
    for i in range(20000):
        acc = f(acc)(i)
    h = hashlib.sha256(repr((sorted(s), top, acc, len(d), sum(d.values()))).encode()).hexdigest()
    return h


def worker(args):
    cpu, seconds = args
    try:
        os.sched_setaffinity(0, {cpu})
        import numpy as np
        from numba import njit

        @njit(cache=False)
        def kern(a, n):
            s = 0.0
            for i in range(n):
                s += np.sin(a[i % a.shape[0]]) * np.exp(-a[(i * 7) % a.shape[0]] * 1e-3)
            return s

        rng = np.random.default_rng(12345)
        a = rng.standard_normal(2048)
        m = rng.standard_normal((384, 384))

        def one():
            h = hashlib.sha256()
            h.update(py_workload(7).encode())
            h.update(str(float(kern(a, 600000))).encode())
            x = m.copy()
            for _ in range(3):
                x = x @ m
                x /= np.abs(x).max()
            h.update(x.tobytes())
            return h.hexdigest()

        ref, reps, bad, errs = None, 0, 0, []
        t0 = time.time()
        while time.time() - t0 < seconds or reps < 2:
            try:
                d = one()
            except Exception as e:      # an interpreter-level anomaly is itself the signal
                errs.append(repr(e)[:120])
                reps += 1
                continue
            if ref is None:
                ref = d
            elif d != ref:
                bad += 1
            reps += 1
        return dict(cpu=cpu, reps=reps, mismatches=bad, exceptions=errs[:3], n_exceptions=len(errs), digest=ref)
    except Exception:
        return dict(cpu=cpu, fatal=traceback.format_exc()[-300:])


if __name__ == '__main__':
    seconds = float(sys.argv[1])
    cpus = [int(c) for c in sys.argv[2].split(',')] if len(sys.argv) > 2 else sorted(os.sched_getaffinity(0))
    os.nice(19)
    with mp.get_context('spawn').Pool(len(cpus)) as pool:
        res = pool.map(worker, [(c, seconds) for c in cpus], chunksize=1)
    digests = Counter(r.get('digest') for r in res)
    print(json.dumps(dict(results=res, distinct_digests=len(digests), digest_counts=digests.most_common()), indent=1))
