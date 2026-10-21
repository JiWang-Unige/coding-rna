"""Sequential queue (one GPU job at a time): for each setting -> m30_run_setting (train + internal validation) -> B inference of the selected checkpoint.
Restartable: finished steps are skipped.  Usage: m30_queue.py <setting> [<setting> ...]   setting = species:k:treatment, e.g. dryas_octopetala:0.25:U"""
import json, os, subprocess, sys, time
from pathlib import Path
import m30_common as M
from m30_run_setting import PY, CPUS, ENV, HERE, sh, gpu_free

logs = Path(os.path.expanduser('~/coding-rna/logs'))
for spec in sys.argv[1:]:
    species, k, treat = spec.split(':')
    name = f'{species}_k{int(round(float(k)*100))}_{treat}'
    out = M.OUT / 'runs' / name
    if not (out / 'selection.json').exists() or 'selected_step' not in json.load(open(out / 'selection.json')):
        rc = subprocess.run([PY, 'm30_run_setting.py', species, k, treat], cwd=HERE, env=ENV, stdin=subprocess.DEVNULL,
                            stdout=open(logs / f'm30-setting-{name}.out', 'ab'), stderr=subprocess.STDOUT).returncode
        if rc != 0:
            print('setting failed', name, rc, flush=True); continue
    sel = json.load(open(out / 'selection.json'))
    step = sel['selected_step']
    evd = M.OUT / 'evalB' / name
    B = ','.join(M.SPECIES[species]['B'])
    for attempt in range(1, 8):
        if list(evd.glob('DONE_*.json')):
            break
        while gpu_free() < 6000:
            time.sleep(120)
        sh([PY, 'm30_infer.py', '--model', str(out / f'model_{step:06d}.pt'), '--seqs', B, '--out', str(evd), '--decode-workers', '4', '--cpus', CPUS],
           logs / f'm30-evalB-{name}.log')
    print('setting complete', name, 'selected step', step, flush=True)
