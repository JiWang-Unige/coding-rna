"""PANEL-P6: run one command in its own process group under a host-RAM guard (stdlib only).  Same interface and result file as M38 `vip_guard.py`, but the memory that is
limited is the summed PSS of the group (/proc/<pid>/smaps_rollup), not the summed RSS: ANNEVO's 8 forked DataLoader workers share their parent's copy-on-write pages, and
an RSS sum counts those pages once per worker (2026-11-09: ANNEVO on Dryas was killed at an RSS sum of 73.7 GB with 119 GB of the host still available).
The RSS sum is still recorded (peak_rss_sum_mb).

  python3.12 -I p6_guard.py --out RES.json [--max-rss-gb 64] [--timeout-s N] [--gpu] [--interval 5] -- COMMAND ARGS...
RES.json: rc, wall_seconds, peak_rss_mb (= peak PSS sum, the guarded quantity), peak_rss_sum_mb, peak_gpu_mb, killed, reason.
"""
import argparse
import json
import os
import signal
import subprocess
import sys
import time


def group_pids(pgid):
    out = subprocess.run(['ps', '-eo', 'pid=,pgid='], capture_output=True, text=True).stdout
    return [int(p[0]) for p in (line.split() for line in out.splitlines()) if len(p) == 2 and int(p[1]) == pgid]


def mem_kb(pid):
    """-> (pss_kb, rss_kb) of one process (0, 0 if it is gone)"""
    pss = rss = 0
    try:
        with open(f'/proc/{pid}/smaps_rollup') as f:
            for line in f:
                if line.startswith('Pss:'):
                    pss = int(line.split()[1])
                elif line.startswith('Rss:'):
                    rss = int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        pass
    return pss, rss


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
    ap.add_argument('--max-rss-gb', type=float, default=64.0, help='limit on the summed PSS of the process group')
    ap.add_argument('--timeout-s', type=float, default=0)
    ap.add_argument('--gpu', action='store_true')
    ap.add_argument('--interval', type=float, default=5.0)
    ap.add_argument('cmd', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd
    t0 = time.time()
    proc = subprocess.Popen(cmd, start_new_session=True)
    pgid = proc.pid
    peak_pss = peak_rss = 0.0
    peak_gpu = 0 if a.gpu else None
    killed, reason, last = False, None, 0.0
    while proc.poll() is None:
        if last == 0.0 or time.time() - last >= a.interval:
            last = time.time()
            pids = group_pids(pgid)
            m = [mem_kb(p) for p in pids]
            pss, rss = sum(x[0] for x in m) / 1024.0, sum(x[1] for x in m) / 1024.0
            peak_pss, peak_rss = max(peak_pss, pss), max(peak_rss, rss)
            if a.gpu and pids:
                g = gpu_mb(set(pids))
                if g is not None:
                    peak_gpu = max(peak_gpu, g)
            if pss > a.max_rss_gb * 1024:
                killed, reason = True, f'PSS {pss / 1024:.1f} GB > {a.max_rss_gb} GB'
        if not killed and a.timeout_s and time.time() - t0 > a.timeout_s:
            killed, reason = True, f'timeout {a.timeout_s}s'
        if killed:
            print(f'[p6_guard] killing process group {pgid}: {reason}', file=sys.stderr, flush=True)
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
    # leftovers of the group (2026-11-10: after a Carbon-A crash its multiprocessing-spawn GPU worker survived, reparented to init, holding 6.7 GB of GPU memory for hours)
    orphans = [p for p in group_pids(pgid) if p != os.getpid()]
    if orphans:
        print(f'[p6_guard] killing {len(orphans)} leftover process(es) of group {pgid}: {orphans}', file=sys.stderr, flush=True)
        try:
            os.killpg(pgid, signal.SIGTERM)
            time.sleep(5)
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    res = dict(cmd=cmd, rc=(137 if killed else rc), wall_seconds=round(time.time() - t0, 1), peak_rss_mb=round(peak_pss, 1), peak_rss_sum_mb=round(peak_rss, 1),
               peak_gpu_mb=peak_gpu, killed=killed, reason=reason, orphans_killed=len(orphans), time=time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(a.out, 'w') as h:
        json.dump(res, h)
    print('[p6_guard] ' + json.dumps({k: v for k, v in res.items() if k != 'cmd'}), file=sys.stderr, flush=True)
    sys.exit(res['rc'])


if __name__ == '__main__':
    main()
