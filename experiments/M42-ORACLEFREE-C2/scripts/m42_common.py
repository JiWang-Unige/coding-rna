"""M42 gate 1 shared setup (protocol 35 §4).  Paths are relative to the repository root (ROOT/data, ROOT/outputs), so the same code runs on ont and on any
machine where HOME holds `coding-rna-reorg` with real `data/` and `outputs/` directories.  No annotation of the target species is read here.
"""
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
ROOT = HERE.parents[2]
SP = os.environ.get('M42_SPECIES', 'ziziphus_jujuba')
DATA = ROOT / 'data' / 'm42' / SP                    # positives.json, background_V*.json, split.json, eval_regions.txt (m42_prep.py)
M1 = ROOT / 'data' / 'm1_screen' / SP                # genome.fa, lengths.json (the path m29_common.read_sequence uses)
OUT = Path(os.environ['M42_OUT']) if os.environ.get('M42_OUT') else ROOT / 'outputs' / 'M42-ORACLEFREE-C2'
BASE = ROOT / 'outputs' / 'M34-NTV2-SCALEUP' / 'pc2_T2' / 'model_002304.pt'   # protocol 35 §4.1: T2 = PlantCAD2 + LoRA, step 2304

# protocol 35 §4.1: treatment of the M36 seed labeler, unknown mask (None = no known background needed), loss
ARMS = {
    'U': dict(treatment='S', mask='V1', loss='crf'),     # partial-observation CRF: positives labelled, V1 unknown, rest background
    'L': dict(treatment='S', mask='V1', loss='ce'),      # same targets, per-position cross-entropy on the known positions only; same CRF decoder
    'F': dict(treatment='F', mask=None, loss='crf'),     # fill background: everything but the positives is background
    'U3': dict(treatment='S', mask='V3', loss='crf'),    # all non-positive unknown except the positives' +-1 kb flanks
}


def split():
    return json.load(open(DATA / 'split.json'))


def lengths():
    return {k: int(v) for k, v in json.load(open(M1 / 'lengths.json')).items()}


def a_lengths():
    L = lengths()
    return {c: L[c] for c in split()['A']}


def eval_regions():
    return (DATA / 'eval_regions.txt').read_text().strip()


def register(MC=None):
    """Make the frozen m29_common reader (read_sequence / LENGTH / SEQ2SPECIES) aware of the target-species chromosomes."""
    if MC is None:
        import m29_common as MC
    for c, n in lengths().items():
        MC.SEQ2SPECIES[c] = SP
        MC.LENGTH[c] = n
    return MC


def load_background(mask):
    raw = json.load(open(DATA / f'background_{mask}.json'))
    return {c: [(int(a), int(b)) for a, b in iv] for c, iv in raw.items()}
