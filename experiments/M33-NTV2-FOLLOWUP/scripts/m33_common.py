"""M33 shared setup: paths and imports (frozen M29/M30/M31 code is imported, never edited)."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[0].parent
for p in (str(EXP / 'M31-C1-NTV2' / 'scripts'), str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np
import m30_common as M30
import a1_common as K
M29 = M30.register()
sys.path.insert(0, str(K.CA_SRC))
OUT = K.OUTS / 'M33-NTV2-FOLLOWUP'
M31OUT = K.OUTS / 'M31-C1-NTV2'
MODEL = OUT / 'model_B_step1536.pt'      # copy of M31 runB/model_001536.pt (kept off command lines: M32's queue driver blocks on any process mentioning the M31 directory name)
HELDOUT = ('NC_003071.7', 'NC_089044.1')          # TRAIN internal-validation chromosomes of M31 (never used to fit B; B checkpoint chosen on them)
DEV = ('NC_003074.8', 'NC_089041.1')
TAU_GRID = (0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)    # declared in protocol.md before any B posterior existed


def lp_by_key(post_dir, tau, seqs):
    """{(seq,strand,exons-tuple): log_p} for all decoded chains of the given contigs at one tau"""
    out = {}
    for s in seqs:
        for st in '+-':
            f = Path(post_dir) / f'tau{tau:g}' / f'{s}_{st}.json'
            for g in json.load(open(f))['genes']:
                out[(s, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
    return out


def calibration(prob, label):
    from m29_r3_posterior_filter import reliability, auroc_ap
    bins, ece = reliability(prob, label)
    auc, ap = auroc_ap(prob, label)
    return dict(n=len(label), hits=int(np.sum(label)), AUROC=auc, AP=ap, ECE=ece, bins=bins)
