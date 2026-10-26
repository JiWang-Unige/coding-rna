"""CPU task queue of the M36 pipeline (pattern of M33 `m33_post_queue.py`: a task starts as soon as the potentials / decoded chains it needs exist, largest contigs first,
bounded by a worker count and a memory budget, pinned to --cpus).  Two modes, both resumable (a task whose output exists is skipped) and single-instance (flock):

  seeds   zero-shot run on the label-set A contigs  ->  per (contig, strand): gene-set posterior (m33_posterior.py, tau = 1) and the per-base 'no gene' marginal
          (m36_nogene_marginal.py); when every task of a species is done, `m36_s1_seeds.py build <species>` writes the S1 seeds / pseudo-background / READY_<species>.
  bpost   evaluation-set-B runs of the trained arms (and Z)  ->  per (contig, strand) posterior only (descriptive posterior-filtered F1, protocol 22 §4);
          when every task of a run is done its pot/ is deleted (never the Z run's).  Nothing here reads any reference annotation.

usage: m36_cpu_queue.py seeds [--run DIR] [--post DIR] [--nogene DIR] [--species a,b] [common options]
       m36_cpu_queue.py bpost --runs NAME=RUNDIR:SPECIES[:keep] [...] [--post-root DIR]            (post-root/<NAME>/tau1/<contig>_<strand>.json)
common: --workers 8 --cpus 24-31 --mem-gb 40 --idle-timeout-h 30 --log FILE --poll 30
"""
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
PY = os.path.expanduser('~/envs/generanno/bin/python')
M33_POST = EXPS / 'M33-NTV2-FOLLOWUP' / 'scripts' / 'm33_posterior.py'
MEM_PER_BASE = dict(post=130, nogene=400)         # bytes per contig base: M33 production ~120 (posterior); m36_nogene_marginal: float32 pot + scaled copy + two float64 [n,15] temporaries (4-5 GB per 20 Mb, 400 B/base kept as the safe side)


class Task:
    def __init__(self, kind, group, run, seq, strand, out, length, cmd, done_file):
        self.kind, self.group, self.run, self.seq, self.strand, self.out = kind, group, Path(run), seq, strand, Path(out)
        self.length, self.cmd, self.done_file = length, cmd, Path(done_file)

    @property
    def name(self):
        return f'{self.kind}:{self.group}:{self.seq}:{self.strand}'

    def done(self):
        return self.done_file.exists()

    def ready(self):
        pot = self.run / 'pot' / f'{self.seq}_{self.strand}.npy'
        if not pot.exists():
            return False
        if self.kind == 'nogene':
            return True
        ch = self.run / 'chains' / f'{self.seq}_{self.strand}.json'
        if not ch.exists():
            return False
        try:
            json.load(open(ch))
        except Exception:
            return False
        return True

    def mem(self):
        return self.length * MEM_PER_BASE[self.kind]


def post_task(group, run, seq, strand, tau, out, length):
    return Task('post', group, run, seq, strand, out, length, [PY, str(M33_POST), str(run), seq, strand, f'{tau:g}', str(out)],
                Path(out) / f'tau{tau:g}' / f'{seq}_{strand}.json')


def nogene_task(group, run, seq, strand, tau, out, length):
    return Task('nogene', group, run, seq, strand, out, length, [PY, str(HERE / 'm36_nogene_marginal.py'), 'run', str(run), seq, strand, f'{tau:g}', str(out)],
                Path(out) / f'tau{tau:g}' / f'nogene_{seq}_{strand}.json')


def run_tasks(tasks, hooks, workers, mem_budget, cpus, emit, idle_timeout=30 * 3600, poll=30., niceness=10, max_tries=2):
    """tasks: [Task]; hooks: {group: callable() | (callable(), already_done_predicate)} fired once when every task of the group is done.  -> (failed task names, hooks run)"""
    pending = [t for t in tasks if not t.done()]
    running, tries, fired = [], {}, set()
    last = time.time()
    emit(event='start', tasks=len(tasks), pending=len(pending))

    def group_done(g):
        return all(t.done() for t in tasks if t.group == g)

    def fire():
        for g, h in hooks.items():
            fn, skip = h if isinstance(h, tuple) else (h, None)
            if g not in fired and group_done(g):
                fired.add(g)
                if skip is not None and skip():
                    emit(event='group_done', group=g, hook='already done')
                    continue
                emit(event='group_done', group=g)
                try:
                    fn()
                except BaseException as e:
                    emit(event='hook_failed', group=g, error=repr(e)[:500])
                    raise
                emit(event='hook_done', group=g)

    fire()
    while pending or running:
        for r in list(running):
            t, p, f = r
            rc = p.poll()
            if rc is None:
                continue
            running.remove(r)
            f.close()
            last = time.time()
            if rc == 0 and t.done():
                emit(event='done', task=t.name)
            else:
                tries[t.name] = tries.get(t.name, 0) + 1
                emit(event='fail', task=t.name, rc=rc, tries=tries[t.name])
                if tries[t.name] < max_tries:
                    pending.append(t)
        fire()
        inflight = sum(t.mem() for t, _, _ in running)
        cand = sorted((t for t in pending if t.ready() and not t.done()), key=lambda t: -t.length)
        for t in cand:
            if len(running) >= workers:
                break
            if running and inflight + t.mem() > mem_budget:
                continue
            log = t.done_file.parent / f'{t.done_file.stem}.out'
            log.parent.mkdir(parents=True, exist_ok=True)
            f = open(log, 'a')
            cmd = ['taskset', '-c', cpus, 'nice', '-n', str(niceness)] + t.cmd if cpus else t.cmd
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
            running.append((t, p, f))
            pending.remove(t)
            inflight += t.mem()
            last = time.time()
            emit(event='launch', task=t.name, bases=t.length, pid=p.pid)
        pending = [t for t in pending if not t.done()]
        if pending and not running and time.time() - last > idle_timeout:
            emit(event='idle_timeout', waiting=len(pending))
            break
        if pending or running:
            time.sleep(poll)
    fire()
    failed = sorted({t.name for t in tasks if not t.done()})
    emit(event='finished', failed=failed)
    return failed, sorted(fired)


def lock(name):
    import m36_common as W
    W.OUT.mkdir(parents=True, exist_ok=True)
    f = open(W.OUT / f'.cpuq_{name}.lock', 'w')
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print(f'cpu queue {name} already running')
        sys.exit(3)
    return f


def contig_lengths(species):
    import m36_common as W
    return W.a_lengths(species)


def seeds_tasks(run, post, nogene, species):
    import m36_common as W
    tasks, hooks = [], {}
    for sp in species:
        for c, L in contig_lengths(sp).items():
            for st in '+-':
                tasks.append(post_task(sp, run, c, st, W.S1_TAU, post, L))
                tasks.append(nogene_task(sp, run, c, st, W.S1_TAU, nogene, L))
    return tasks


def seeds_hook(species, post, nogene, seeds_out, cpus):
    def fn():
        cmd = ['taskset', '-c', cpus, 'nice', '-n', '10', PY, str(HERE / 'm36_s1_seeds.py'), 'build', species, '--post', str(post), '--nogene', str(nogene), '--out', str(seeds_out)]
        r = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True)
        (Path(seeds_out) / f'build_{species}.log').write_text(r.stdout + '\n' + r.stderr)
        if r.returncode != 0:
            raise RuntimeError(f'm36_s1_seeds build {species} failed rc={r.returncode}: {r.stderr[-800:]}')
    return fn


def bpost_tasks(runs, post_root):
    """runs: [(name, run_dir, species, keep_pot)]"""
    import m36_common as W
    import m30_common as M30
    tasks, hooks = [], {}
    lens_all = {}
    for name, run, sp, keep in runs:
        lens = lens_all.setdefault(sp, M30.lengths(sp))
        for c in M30.SPECIES[sp]['B']:
            for st in '+-':
                tasks.append(post_task(name, run, c, st, W.S1_TAU, Path(post_root) / name, lens[c]))
        if not keep:
            hooks[name] = (lambda r=Path(run): shutil.rmtree(r / 'pot', ignore_errors=True))
    return tasks, hooks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('seeds', 'bpost'))
    ap.add_argument('--run')
    ap.add_argument('--post')
    ap.add_argument('--nogene')
    ap.add_argument('--seeds-out')
    ap.add_argument('--species', default=None)
    ap.add_argument('--runs', nargs='*', default=[])
    ap.add_argument('--post-root')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--cpus', default='24-31')
    ap.add_argument('--mem-gb', type=float, default=40.)
    ap.add_argument('--idle-timeout-h', type=float, default=30.)
    ap.add_argument('--poll', type=float, default=30.)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    import m36_common as W
    W.K.repoint()
    lk = lock(a.mode)
    log = open(a.log, 'a') if a.log else sys.stdout
    emit = lambda **k: (log.write(json.dumps(dict(k, t=time.strftime('%F %T'))) + '\n'), log.flush())
    if a.mode == 'seeds':
        species = a.species.split(',') if a.species else list(W.SPECIES)
        run, post, nogene, out = (Path(a.run or W.OUT / 'zeroshot_A'), Path(a.post or W.OUT / 'post_A'), Path(a.nogene or W.OUT / 'nogene_A'),
                                  Path(a.seeds_out or W.OUT / 'seeds'))
        tasks = seeds_tasks(run, post, nogene, species)
        hooks = {sp: (seeds_hook(sp, post, nogene, out, a.cpus), (lambda sp=sp: (out / f'READY_{sp}').exists())) for sp in species}
    else:
        runs = []
        for s in a.runs:
            name, _, rest = s.partition('=')
            parts = rest.split(':')
            runs.append((name, parts[0], parts[1], len(parts) > 2 and parts[2] == 'keep'))
        tasks, hooks = bpost_tasks(runs, a.post_root or W.OUT / 'postB')
    failed, fired = run_tasks(tasks, hooks, a.workers, a.mem_gb * 2**30, a.cpus, emit, a.idle_timeout_h * 3600, a.poll)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
