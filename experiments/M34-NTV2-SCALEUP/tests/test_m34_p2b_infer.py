"""P2b inference wrapper (m34_p2b_infer.py), CPU only: the M30 evaluation-set-B contigs must be readable through the registration used by m34_infer (main process
and spawned decode worker), next to the M34 validation pseudo-contigs.
Run on ont:  cd .../M34-NTV2-SCALEUP/tests && M34_BACKEND=pc2 taskset -c 25-31 nice -n 10 ~/envs/plantcad2/bin/python -m pytest -q test_m34_p2b_infer.py
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import m34_train_common as T

pytestmark = pytest.mark.skipif(T.BACKEND != 'pc2', reason='run with M34_BACKEND=pc2 under ~/envs/plantcad2')

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
B = {'dryas_octopetala': ['NC_141367.1', 'NC_141369.1', 'NC_141371.1', 'NC_141373.1'],
     'eschscholzia_californica': ['NW_028238309.1', 'NW_028238461.1', 'NW_028238483.1']}
ALL_B = [c for cs in B.values() for c in cs]


def run_py(code):
    return subprocess.run([sys.executable, '-c', code], cwd=SCRIPTS, capture_output=True, text=True, env=os.environ.copy())


def test_unwrapped_registration_does_not_know_the_B_contigs():
    r = run_py("import m34_train_common as T; MC = T.register(); MC.LENGTH['NC_141367.1']")
    assert r.returncode != 0 and 'KeyError' in r.stderr


def test_wrapper_registers_B_contigs_and_keeps_m34_ones():
    import m34_p2b_infer as W
    import m30_common as M30
    assert W.X.T.register is W._register_all and T.register is W._register_all
    assert sys.path[0] == str(SCRIPTS) and sys.path[1] == str(T.ENC_SCRIPTS)            # the backend's script dirs stay first
    MC = T.register()
    for sp, cs in B.items():
        assert sorted(M30.SPECIES[sp]['B']) == sorted(cs)                              # protocol 13 section 2 evaluation set B
        L = M30.lengths(sp)
        for c in cs:
            assert MC.LENGTH[c] == L[c] > 1_000_000 and MC.SEQ2SPECIES[c] == sp
    k = T.new_keys()[0]
    v = T.val_contig(k)
    assert MC.LENGTH[v] == 8_000_000 and MC.SEQ2SPECIES[v] == k                         # M34 pseudo-contigs still registered
    assert len(MC.read_sequence(v)) == 8_000_000
    for c in ('NC_141373.1', 'NW_028238483.1'):                                         # one contig per species is read in full (the others share the same code path)
        s = MC.read_sequence(c)
        assert len(s) == MC.LENGTH[c] and set(s) <= set('ACGTN')
        assert W.X.I.resolve(c, 0, None, len(s)) == (len(s), c)
    assert W.X.I.parse_region('NC_141367.1') == ('NC_141367.1', 0, None)


def test_spawned_decode_worker_sees_the_B_contigs():
    r = subprocess.run([sys.executable, 'm34_p2b_infer.py', '--probe', 'NC_141373.1,NW_028238483.1'], cwd=SCRIPTS, capture_output=True, text=True, env=os.environ.copy())
    assert r.returncode == 0, r.stderr[-2000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out['backend'] == 'pc2' and out['same'] is True
    assert [x['contig'] for x in out['worker']] == ['NC_141373.1', 'NW_028238483.1']
    assert all(x['seq_len'] == x['length'] for x in out['worker'])
