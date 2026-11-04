"""M38-VIPSANIA: append one JSON event line (or write a JSON file) from key=value arguments (stdlib only).
  vip_event.py --file events.jsonl event=step_done step=ara_pre [--merge guard_result.json] [--write] key=value ...
Values that look like int / float / true / false / null are typed.  --merge copies rc, wall_seconds (as `seconds`), peak_rss_mb, peak_gpu_mb, killed, reason from a vip_guard.py
result file.  --write overwrites FILE with the single object (pretty-printed) instead of appending a line."""
import argparse
import json
import time


def typed(v):
    if v in ('true', 'false'):
        return v == 'true'
    if v == 'null':
        return None
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--file', required=True)
    ap.add_argument('--merge', default=None)
    ap.add_argument('--write', action='store_true')
    ap.add_argument('kv', nargs='*')
    a = ap.parse_args()
    ev = {}
    for item in a.kv:
        k, _, v = item.partition('=')
        ev[k] = typed(v)
    ev['time'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    if a.merge:
        try:
            g = json.load(open(a.merge))
            ev.update(rc=g.get('rc'), seconds=g.get('wall_seconds'), peak_rss_mb=g.get('peak_rss_mb'), peak_gpu_mb=g.get('peak_gpu_mb'), killed=g.get('killed'))
            if g.get('reason'):
                ev['reason'] = g['reason']
        except Exception as e:
            ev['merge_error'] = str(e)
    if a.write:
        json.dump(ev, open(a.file, 'w'), indent=1)
    else:
        with open(a.file, 'a') as h:
            h.write(json.dumps(ev) + '\n')


if __name__ == '__main__':
    main()
