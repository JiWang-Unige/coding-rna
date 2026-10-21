"""Dry run of m33_select.py / m33_dev.py / m33_evalB.py on EXISTING outputs with FAKE posteriors (monotone transform of the decode gain) to exercise the code paths
(file layout, key matching, scoring, bootstrap).  Numbers are meaningless; the evalB part uses the M30 gate chains and posteriors as stand-in for model B, so its
NTv2B_post arm must reproduce the Gate_post arm exactly (74.09 / 72.52).  Writes only under outputs/M33-NTV2-FOLLOWUP/_dryrun/.
Usage: dryrun_fake_posteriors.py"""
import json, subprocess, sys, os
from pathlib import Path
HERE = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(HERE))
import numpy as np
import m33_common as C
D = C.OUT / '_dryrun'
PY = os.path.expanduser('~/envs/generanno/bin/python')
M31 = next(C.K.OUTS.glob('M31*'))
def fake(run, seqs, taus, out):
    for t in taus:
        for s in seqs:
            for st in '+-':
                d = json.load(open(Path(run) / 'chains' / f'{s}_{st}.json'))
                rows = [dict(exons=c['exons'], gain=c['gain'], log_p=float(min(0.0, -2 * np.exp(-t * c['gain'] / 10)))) for c in d['chains']]
                f = Path(out) / f'tau{t:g}' / f'{s}_{st}.json'
                f.parent.mkdir(parents=True, exist_ok=True)
                json.dump(dict(seq=s, strand=st, tau=t, genes=rows, log_partition=0.0, seconds=0), open(f, 'w'))
fake(M31 / 'runB' / 'val_1536', C.HELDOUT, (0.5, 1.0), D / 'post_val')
fake(M31 / 'devB_final', C.DEV, (0.5, 1.0), D / 'post_dev')
def run(*a):
    print('$', ' '.join(map(str, a)), flush=True)
    subprocess.run([PY, *map(str, a)], check=True, cwd=HERE)
run(HERE / 'm33_select.py', M31 / 'runB' / 'val_1536', D / 'post_val', D / 'selection.json')
run(HERE / 'm33_dev.py', M31 / 'devB_final', D / 'post_dev', D / 'selection.json', D / 'dev_report.json')
# evalB with the gate as stand-in: selection tau=0.03 / thr = M30 threshold, posterior dir = M30 gate B posteriors
import m30_common as M30
sel = dict(selection=dict(tau=0.03, thr=json.load(open(M30.OUT / 'gate_posterior_threshold.json'))['threshold']))
json.dump(sel, open(D / 'selection_gate.json', 'w'))
pd = D / 'post_gate_B'
for t in (pd / 'tau0.03',):
    t.parent.mkdir(parents=True, exist_ok=True)
    if not t.exists():
        os.symlink(M30.OUT / 'posterior_gate_tau0.03_B', t)
run(HERE / 'm33_evalB.py', 'dryas_octopetala', M30.OUT / 'evalB' / 'gate', pd, D / 'selection_gate.json', D / 'evalB_dryas_fake.json')
