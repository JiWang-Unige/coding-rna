"""M38-VIPSANIA: run one command in its own process group with a host-RAM guard and resource accounting (stdlib only).

  python3.12 -I vip_guard.py --out RES.json [--max-rss-gb 40] [--timeout-s N] [--gpu] [--interval 5] -- COMMAND ARGS...

Polls the summed RSS of the whole process group (and, with --gpu, the GPU memory that nvidia-smi attributes to its PIDs).  If the RSS exceeds --max-rss-gb (or the
timeout expires) the group is killed (SIGTERM, then SIGKILL after 20 s), so that a runaway CPU job can never push the shared server into the OOM killer.  The exit code of
the command is the exit code of this wrapper (137 if the guard killed it).  RES.json: rc, wall_seconds, peak_rss_mb, peak_gpu_mb (null without --gpu), killed, reason.
stdout/stderr/stdin of the command are inherited unchanged.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import time


def group_procs(pgid):
    out = subprocess.run(['ps', '-eo', 'pid=,pgid=,rss='], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        p = line.split()
        if len(p) == 3 and int(p[1]) == pgid:
            rows.append((int(p[0]), int(p[2])))
    return rows


def gpu_mb(pids):
    try:
        out = subprocess.run(['nvidia-smi', '--query-compute-apps=pid,used_memory', '--format=csv,noheader,nounits'], capture_output=True, text=True, timeout=20).stdout
    except Exception:
        return None
    tot = 0
    for line in out.splitlines():
        p = [x.strip() for x in line.split(',')]
        if len(p) == 2 and p[0].isdigit() and int(p[0]) in pids and p[1].isdigit():
            tot += int(p[1])
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-rss-gb', type=float, default=40.0)
    ap.add_argument('--timeout-s', type=float, default=0)
    ap.add_argument('--gpu', action='store_true')
    ap.add_argument('--interval', type=float, default=5.0)
    ap.add_argument('cmd', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd
    t0 = time.time()
    proc = subprocess.Popen(cmd, start_new_session=True)
    pgid = proc.pid
    peak_rss, peak_gpu, killed, reason = 0, (0 if a.gpu else None), False, None
    last = 0.0
    while proc.poll() is None:
        if time.time() - last >= a.interval or last == 0.0:
            last = time.time()
            rows = group_procs(pgid)
            rss = sum(r for _, r in rows) / 1024.0
            peak_rss = max(peak_rss, rss)
            if a.gpu and rows:
                g = gpu_mb({p for p, _ in rows})
                if g is not None:
                    peak_gpu = max(peak_gpu, g)
            if rss > a.max_rss_gb * 1024:
                killed, reason = True, f'RSS {rss / 1024:.1f} GB > {a.max_rss_gb} GB'
        if not killed and a.timeout_s and time.time() - t0 > a.timeout_s:
            killed, reason = True, f'timeout {a.timeout_s}s'
        if killed:
            print(f'[vip_guard] killing process group {pgid}: {reason}', file=sys.stderr, flush=True)
            try:
                os.killpg(pgid, signal.SIGTERM)
                for _ in range(40):
                    if proc.poll() is not None:
                        break
                    time.sleep(0.5)
                os.killpg(pgid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            break
        time.sleep(0.5)
    rc = proc.wait()
    res = dict(cmd=cmd, rc=(137 if killed else rc), wall_seconds=round(time.time() - t0, 1), peak_rss_mb=round(peak_rss, 1),
               peak_gpu_mb=peak_gpu, killed=killed, reason=reason, time=time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(a.out, 'w') as h:
        json.dump(res, h)
    print('[vip_guard] ' + json.dumps({k: v for k, v in res.items() if k != 'cmd'}), file=sys.stderr, flush=True)
    sys.exit(res['rc'])


if __name__ == '__main__':
    main()
