"""M29-A1 shared setup.  Frozen code (pipeline/m28, refs/carbon-a) is only imported, never edited.
Run:  cd ~/coding-rna-reorg && taskset -c <cpus> ~/envs/generanno/bin/python experiments/M29-A1-RESIDUAL/scripts/<script>.py
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                                   # ~/coding-rna-reorg
OUTS = ROOT / 'outputs'                                  # -> ~/coding-rna/outputs
A1 = OUTS / 'M29-A1-RESIDUAL'
M29 = OUTS / 'M29-CARBON-A-GAP'
CACHE = M29 / '_cache'
CA_SRC = ROOT / 'experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/src'
EXTRA = os.path.expanduser('~/envs/m29_extra')
for p in (str(ROOT / 'pipeline/m28'), str(ROOT / 'pipeline/m28/scripts'),
          str(ROOT / 'experiments/M29-CARBON-A-GAP/scripts'), str(HERE), EXTRA):
    if p not in sys.path:
        sys.path.insert(0, p)
BAOBAB = '/home/users/j/jwang/coding-rna'


def repoint():
    """Some frozen modules hard-code the Baobab ROOT; re-point module-level Paths (no frozen file is edited)."""
    for m in list(sys.modules.values()):
        for k, v in list(getattr(m, '__dict__', {}).items()):
            if isinstance(v, Path) and str(v).startswith(BAOBAB):
                setattr(m, k, ROOT / Path(str(v)).relative_to(BAOBAB))


def pin(cpus):
    os.sched_setaffinity(0, set(cpus))
