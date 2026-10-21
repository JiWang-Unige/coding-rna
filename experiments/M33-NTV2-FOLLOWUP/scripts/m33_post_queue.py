"""CPU queue for `m33_posterior.py`: one subprocess per (run, seq, strand, tau), started as soon as the potentials (pot/*.npy, atomically renamed) and the
decoded chains (chains/*.json, must parse) of that (seq, strand) exist; largest contigs first; at most --workers concurrent processes and at most --max-bases
(sum of contig lengths of running tasks, memory guard ~120 B/base).  Pins to --cpus (never CPU 0/1).  Exits when everything is done or --idle-timeout seconds pass
without progress while tasks are still waiting for their inputs.
Usage: m33_post_queue.py --spec RUN:SEQ,SEQ:TAU,TAU:OUTDIR [--spec ...] --workers 8 --cpus 2,4,... --log FILE
"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
PY = os.path.expanduser('~/envs/generanno/bin/python')


def ready(run, seq, strand):
    pot = Path(run) / 'pot' / f'{seq}_{strand}.npy'
    ch = Path(run) / 'chains' / f'{seq}_{strand}.json'
    if not (pot.exists() and ch.exists()):
        return False
    try:
        json.load(open(ch))
    except Exception:
        return False
    return True


def length(run, seq, strand):
    import numpy as np
    return int(np.load(Path(run) / 'pot' / f'{seq}_{strand}.npy', mmap_mode='r').shape[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spec', action='append', required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--cpus', default='2,4,6,8,10,12,14,16,18,20,22,24')
    ap.add_argument('--max-bases', type=float, default=400e6)
    ap.add_argument('--idle-timeout', type=float, default=6 * 3600)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    tasks = []
    for s in a.spec:
        run, seqs, taus, outd = s.split(':')
        for seq in seqs.split(','):
            for st in '+-':
                for tau in taus.split(','):
                    tasks.append(dict(run=run, seq=seq, strand=st, tau=tau, out=outd))
    done = lambda t: (Path(t['out']) / f"tau{float(t['tau']):g}" / f"{t['seq']}_{t['strand']}.json").exists()
    pending = [t for t in tasks if not done(t)]
    running = []                                   # (task, Popen, bases)
    fails = {}
    last = time.time()
    log = open(a.log, 'a') if a.log else sys.stdout
    emit = lambda **k: (log.write(json.dumps(dict(k, t=time.strftime('%F %T'))) + '\n'), log.flush())
    cpus = a.cpus
    emit(event='start', pending=len(pending))
    while pending or running:
        for r in list(running):
            rc = r[1].poll()
            if rc is None:
                continue
            running.remove(r)
            last = time.time()
            if rc == 0 and done(r[0]):
                emit(event='done', **{k: r[0][k] for k in ('seq', 'strand', 'tau')})
            else:
                key = (r[0]['run'], r[0]['seq'], r[0]['strand'], r[0]['tau'])
                fails[key] = fails.get(key, 0) + 1
                emit(event='fail', rc=rc, tries=fails[key], **{k: r[0][k] for k in ('seq', 'strand', 'tau')})
                if fails[key] < 2:
                    pending.append(r[0])
        inflight = sum(r[2] for r in running)
        cand = []
        for t in pending:
            if ready(t['run'], t['seq'], t['strand']):
                cand.append((length(t['run'], t['seq'], t['strand']), t))
        cand.sort(key=lambda x: -x[0])
        for L, t in cand:
            if len(running) >= a.workers:
                break
            if inflight + L > a.max_bases and running:
                continue
            p = subprocess.Popen(['taskset', '-c', cpus, PY, str(HERE / 'm33_posterior.py'), t['run'], t['seq'], t['strand'], t['tau'], t['out']],
                                 stdout=open(os.path.join(os.path.dirname(a.log or '/tmp/x'), f"post_{t['seq']}_{t['strand']}_{t['tau']}.out"), 'a'), stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
            running.append((t, p, L)); inflight += L; pending.remove(t); last = time.time()
            emit(event='launch', seq=t['seq'], strand=t['strand'], tau=t['tau'], bases=L)
        if pending and not running and time.time() - last > a.idle_timeout:
            emit(event='idle_timeout', waiting=len(pending))
            break
        time.sleep(30)
    emit(event='finished', failed=[list(k) for k, v in fails.items() if v >= 2])


if __name__ == '__main__':
    main()
