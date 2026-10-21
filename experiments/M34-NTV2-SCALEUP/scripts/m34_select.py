"""Checkpoint selection / one-shot DEV sentinel of an M34 arm.
  m34_select.py select <arm>      pick the final checkpoint = argmax over steps {1536, 2304, 3072} of the Ath/Osa pooled internal-validation F1 (ties -> larger step);
                                  uses only metrics.jsonl lines with tag 'select'; writes best_step.txt + selection.json (refuses to overwrite)
  m34_select.py dev_start <arm>   DEV sentinel: first call records the model (path, size, mtime, md5) in dev_started.json; later calls only succeed for the identical model (inference
                                  itself is resumable, a different model or a finished evaluation is refused).  Exit 0 = may run DEV inference, 11 = refused, 10 = already evaluated.
"""
import hashlib
import json
import sys
import time

import m34_train_common as T

CANDIDATES = (1536, 2304, 3072)


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def main():
    cmd, arm = sys.argv[1], sys.argv[2]
    O = T.arm_dir(arm)
    if cmd == 'select':
        sel = O / 'selection.json'
        if sel.exists():
            print(open(O / 'best_step.txt').read().strip())
            return 0
        rows = [json.loads(l) for l in open(O / 'metrics.jsonl') if l.strip()]
        f1 = {r['step']: r['pooled_F1'] for r in rows if r.get('tag') == 'select' and r.get('arm') == arm}
        cand = {s: f1[s] for s in CANDIDATES if s in f1}
        assert len(cand) == len(CANDIDATES), (cand, 'missing validations')
        best = max(CANDIDATES, key=lambda s: (round(cand[s], 9), s))
        sel.write_text(json.dumps(dict(arm=arm, rule='argmax pooled F1 of Ath+Osa held-out chromosomes over {1536,2304,3072}; ties -> larger step',
                                       candidates=cand, all_validations_select=f1, best_step=best, time=time.strftime('%Y-%m-%dT%H:%M:%S')), indent=1))
        (O / 'best_step.txt').write_text(f'{best}\n')
        print(best)
        return 0
    if cmd == 'dev_start':
        if (O / 'DEV_DONE').exists():
            print('DEV already evaluated')
            return 10
        best = int((O / 'best_step.txt').read_text())
        model = O / f'model_{best:06d}.pt'
        info = dict(arm=arm, step=best, model=str(model), size=model.stat().st_size, md5=md5(model))
        p = O / 'dev_started.json'
        if p.exists():
            old = json.loads(p.read_text())
            if (old['md5'], old['step']) != (info['md5'], info['step']):
                print('DEV was started with a different model', old, info)
                return 11
            return 0
        p.write_text(json.dumps(dict(info, time=time.strftime('%Y-%m-%dT%H:%M:%S')), indent=1))
        return 0
    raise SystemExit('unknown command')


if __name__ == '__main__':
    sys.exit(main())
