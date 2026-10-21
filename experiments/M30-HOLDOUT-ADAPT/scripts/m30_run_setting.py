"""One adaptation setting end to end (protocol §4): train in 256-step stages, validate each checkpoint on the internal validation contig
(labeled-F1), patience 2, select the best trained checkpoint.  Crash-resilient (resume from ckpt_latest, interruptions logged).
Usage: m30_run_setting.py <species> <k> <U|F|FULL>"""
import json, os, subprocess, sys, time
from pathlib import Path
import m30_common as M

CPUS = '2,4,6,8,10,12,14,16-24'
PY = os.path.expanduser('~/envs/generanno/bin/python')
HERE = Path(__file__).resolve().parent
ENV = dict(os.environ, HF_HUB_OFFLINE='1', TRANSFORMERS_OFFLINE='1', TOKENIZERS_PARALLELISM='false', OMP_NUM_THREADS='2', PYTHONUNBUFFERED='1')


def log(out, **kw):
    with open(out / 'interruptions.jsonl', 'a') as f:
        f.write(json.dumps(dict(kw, time=time.strftime('%Y-%m-%dT%H:%M:%S%z'))) + '\n')


def gpu_free():
    return int(subprocess.run(['nvidia-smi', '--query-gpu=memory.free', '--format=csv,noheader,nounits'], capture_output=True, text=True).stdout.split()[0])


def sh(cmd, out_log):
    with open(out_log, 'ab') as o:
        return subprocess.run(['taskset', '-c', CPUS] + cmd, cwd=HERE, env=ENV, stdout=o, stderr=o, stdin=subprocess.DEVNULL).returncode


def main():
    species, k, treat = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    name = f'{species}_k{int(round(k*100))}_{treat}'
    out = M.OUT / 'runs' / name; out.mkdir(parents=True, exist_ok=True)
    logs = Path(os.path.expanduser('~/coding-rna/logs'))
    sel_path = out / 'selection.json'
    sel = json.load(open(sel_path)) if sel_path.exists() else dict(history=[], best=None, stopped=None)
    contig = M.SPECIES[species]['val_contig']
    done_steps = {h['step'] for h in sel['history']}
    for stop in (256, 512, 768, 1024):
        if stop in done_steps:
            continue
        for attempt in range(1, 8):
            while gpu_free() < 9000:
                log(out, event='gpu_busy_wait'); time.sleep(120)
            log(out, event='train_launch', stop=stop, attempt=attempt)
            rc = sh([PY, 'm30_train.py', '--species', species, '--k', str(k), '--treatment', treat, '--out', str(out), '--steps', '1024', '--stop-at', str(stop),
                     '--workers', '12', '--cpus', CPUS, '--resume'], logs / f'm30-train-{name}.log')
            if (out / f'DONE_{stop}.json').exists():
                break
            log(out, event='train_interruption', stop=stop, attempt=attempt, rc=rc)
            time.sleep(30)
        else:
            raise SystemExit('training failed repeatedly')
        vdir = out / f'val_{stop}'
        for attempt in range(1, 8):
            while gpu_free() < 6000:
                time.sleep(120)
            rc = sh([PY, 'm30_infer.py', '--model', str(out / f'model_{stop:06d}.pt'), '--seqs', contig, '--out', str(vdir), '--decode-workers', '4', '--cpus', CPUS],
                    logs / f'm30-valinfer-{name}.log')
            if list(vdir.glob('DONE_*.json')):
                break
            log(out, event='infer_interruption', stop=stop, attempt=attempt, rc=rc); time.sleep(30)
        else:
            raise SystemExit('validation inference failed repeatedly')
        sh([PY, 'm30_valscore.py', str(vdir), species, str(k)], logs / f'm30-valscore-{name}.log')
        v = json.load(open(vdir / 'valscore.json'))
        sel['history'].append(dict(step=stop, labeled_F1=v['F1'], tp=v['tp'], pred=v['pred'], ref=v['ref']))
        best = max(sel['history'], key=lambda h: (h['labeled_F1'], -h['step']))
        sel['best'] = best
        # patience 2: stop after two consecutive checkpoints that do not beat the best so far
        hist = sel['history']
        since = len(hist) - 1 - hist.index(best)
        json.dump(sel, open(sel_path, 'w'), indent=1)
        log(out, event='validated', stop=stop, labeled_F1=v['F1'], best_step=best['step'])
        if since >= 2:
            sel['stopped'] = f'patience2_after_{stop}'; break
    sel['selected_step'] = sel['best']['step']
    json.dump(sel, open(sel_path, 'w'), indent=1)
    print(json.dumps(sel))


if __name__ == '__main__':
    main()
