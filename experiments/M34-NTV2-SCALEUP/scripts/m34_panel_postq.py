"""P5 panel posterior CPU queue: the M36 queue engine (`m36_cpu_queue.Task` / `run_tasks`, imported, not edited) over the panel inference runs.
Per (species, contig, strand): gene-set posterior, tau = 1 (m34_panel_post.py) as soon as the run's pot/<seq>_<strand>.npy and chains/<seq>_<strand>.json exist,
largest contigs first, bounded by --workers and --mem-gb.  When every task of a species is done, that species' pot/ is deleted (the potentials take ~57 bytes
per base and strand, ~0.55 TB for the core 10 species, so they are not kept).  Resumable (a task whose output exists is skipped), single instance (flock).
Nothing here reads a reference annotation.
usage: m34_panel_postq.py --run-root RUNS --post-root POST [--species a,b] [--workers 4 --cpus 24-31 --mem-gb 30 --idle-timeout-h 48 --log FILE]
       runs: RUNS/<species>/{pot,chains}   output: POST/<species>/tau1/<seq>_<strand>.json
"""
import argparse
import fcntl
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M36-BOOTSTRAP-ADAPT' / 'scripts'))
import m36_cpu_queue as Q                                # noqa: E402
import m34_panel_common as PC                            # noqa: E402

TAU = 1


def tasks_for(run_root, post_root, species):
    tasks, hooks = [], {}
    for key in species:
        run, out = Path(run_root) / key, Path(post_root) / key
        L = PC.lengths(key)
        for c in PC.contigs(key):
            for st in '+-':
                tasks.append(Q.Task('post', key, run, c, st, out, L[c], [Q.PY, str(HERE / 'm34_panel_post.py'), str(run), c, st, f'{TAU:g}', str(out)],
                                    out / f'tau{TAU:g}' / f'{c}_{st}.json'))
        hooks[key] = (lambda r=run: shutil.rmtree(r / 'pot', ignore_errors=True))
    return tasks, hooks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-root', required=True)
    ap.add_argument('--post-root', required=True)
    ap.add_argument('--species', default=','.join(PC.CORE10))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--cpus', default='24-31')
    ap.add_argument('--mem-gb', type=float, default=30.)
    ap.add_argument('--idle-timeout-h', type=float, default=48.)
    ap.add_argument('--poll', type=float, default=60.)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    Path(a.post_root).mkdir(parents=True, exist_ok=True)
    lk = open(Path(a.post_root) / '.postq.lock', 'w')
    try:
        fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print('panel posterior queue already running')
        sys.exit(3)
    log = open(a.log, 'a') if a.log else sys.stdout
    emit = lambda **k: (log.write(json.dumps(dict(k, t=time.strftime('%F %T'))) + '\n'), log.flush())
    tasks, hooks = tasks_for(a.run_root, a.post_root, a.species.split(','))
    failed, fired = Q.run_tasks(tasks, hooks, a.workers, a.mem_gb * 2**30, a.cpus, emit, a.idle_timeout_h * 3600, a.poll)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
