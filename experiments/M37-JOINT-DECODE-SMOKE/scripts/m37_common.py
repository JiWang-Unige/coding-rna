"""M37-JOINT-DECODE-SMOKE shared setup (protocol: docs/review_20261016/23_M37_joint_decoding_smoke_protocol.md).
Frozen M29/M30/M31/M33 code is imported, never edited.  Assembly A = Nipponbare AGIS-1.0 NC_089044.1 (chromosome 10; M31 internal
validation, M33 potentials reused).  Assembly B = 93-11 (indica) GCA_000004655.2 ASM465v1 chromosome 10 = CM000135.2 (first candidate of
protocol §1 with an NCBI annotation: MH63RS3 / ZS97RS3 / R498 have none; 9311 ASM465v1 has 37,358 protein-coding genes).
"""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]                                     # .../experiments
for p in (str(EXP / 'M33-NTV2-FOLLOWUP' / 'scripts'), str(EXP / 'M31-C1-NTV2' / 'scripts'), str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np
import m30_common as M30
import a1_common as K
M29 = M30.register()
if str(K.CA_SRC) not in sys.path:
    sys.path.insert(0, str(K.CA_SRC))
OUT = K.OUTS / 'M37-JOINT-DECODE-SMOKE'
DATA = K.ROOT / 'data' / 'm37'                            # -> ~/coding-rna/data/m37
M33OUT = K.OUTS / 'M33-NTV2-FOLLOWUP'
MODEL = K.OUTS / 'M31-C1-NTV2' / 'runB' / 'model_001536.pt'   # frozen, md5 4a89a856...
A_CONTIG = 'NC_089044.1'
B_CONTIG = 'CM000135.2'
B_SPECIES = 'oryza_sativa_9311'                           # data/m37/oryza_sativa_9311/{genome.fa, reference.gff3, lengths.json}
B_ASM = 'GCA_000004655.2_ASM465v1'
B_SPECIES_REL = str(Path('..') / 'm37' / B_SPECIES)       # relative to the frozen M29 DATA root (data/m1_screen) -> no frozen file edited
A_POT = {st: M33OUT / 'val' / 'pot' / f'{A_CONTIG}_{st}.npy' for st in '+-'}
A_CHAINS = {st: M33OUT / 'val' / 'chains' / f'{A_CONTIG}_{st}.json' for st in '+-'}
LAMBDAS = (0.25, 0.5)                                     # protocol §4.2
BLK = 1_000_000                                           # 1-Mb block bootstrap (M31 c1_boot / M33 m33_dev)
MIN_IDENTITY, MIN_BLOCK = 0.90, 1000                      # protocol §3


def register():
    """Make the frozen M29 helpers (sequence reader, reference builder) aware of the B contig; reference pickles go to the M37 cache."""
    p = DATA / B_SPECIES / 'lengths.json'
    if p.exists():
        for k, v in json.load(open(p)).items():
            M29.SEQ2SPECIES[k] = B_SPECIES_REL
            M29.LENGTH[k] = int(v)
    M29.CACHE = OUT / '_cache'
    return M29


def load_chain_json(path):
    """m33_infer chains json -> list of (exons tuple((a,b),...), gain); exons genomic 0-based half-open, ascending"""
    d = json.load(open(path))
    return [(tuple((int(a), int(b)) for a, b in c['exons']), float(c['gain'])) for c in d['chains']], d


def chain_dict(seq, chains_by_strand, cds96=True):
    """{(seq, strand, exons): dict(nex, exons, strand, gain)} unique exact chains, CDS>=96 rule as the project ruler"""
    out = {}
    for st, chains in chains_by_strand.items():
        for ex, gain in chains:
            if cds96 and sum(b - a for a, b in ex) < 96:
                continue
            out.setdefault((seq, st, ex), dict(nex=len(ex), exons=[tuple(e) for e in ex], strand=st, gain=gain))
    return out


def ref_keys(seq):
    M = register()
    ref = M.build_reference(seq)
    keys = {M.key(t) for t in ref['refs']}
    nex = {M.key(t): len(t['CDS']) for t in ref['refs']}
    return ref, keys, nex


def f1(tp, p, r):
    return 100 * 2 * tp / (p + r) if p + r else 0.


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'
