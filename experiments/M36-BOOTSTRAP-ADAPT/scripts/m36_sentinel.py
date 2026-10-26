"""One-shot sentinels of the M36 evaluation-set-B steps (protocol 22 §8: B is evaluated once).  Pure stdlib (no torch / no reference): safe to call from shell drivers.
  start <sentinel.json> <identity-json>   first call atomically creates the sentinel with the identity (what is being run: model md5s, inputs) -> exit 0.
                                          Later calls succeed (exit 0) only for the IDENTICAL identity (the step itself is resumable: a crash may be restarted with the same
                                          model / inputs, never with another one) -> else exit 11.   If <sentinel>.done exists the step already finished -> exit 10.
  done  <sentinel.json>                   mark the step finished (writes <sentinel>.done).
Exit codes:  0 may run   10 already finished (do not run again)   11 refused (different identity than the first start)
"""
import json
import os
import sys
import time


def _canon(x):
    return json.loads(json.dumps(x, sort_keys=True))


def start(path, identity):
    path = str(path)
    if os.path.exists(path + '.done'):
        return 10
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    rec = dict(identity=_canon(identity), time=time.strftime('%Y-%m-%dT%H:%M:%S'), pid=os.getpid())
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        old = json.load(open(path))
        return 0 if old['identity'] == rec['identity'] else 11
    with os.fdopen(fd, 'w') as f:
        json.dump(rec, f, indent=1)
    return 0


def done(path):
    path = str(path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path + '.done', 'w') as f:
        f.write(time.strftime('%Y-%m-%dT%H:%M:%S'))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'start':
        sys.exit(start(sys.argv[2], json.loads(sys.argv[3])))
    if cmd == 'done':
        done(sys.argv[2])
        sys.exit(0)
    raise SystemExit('unknown command')
