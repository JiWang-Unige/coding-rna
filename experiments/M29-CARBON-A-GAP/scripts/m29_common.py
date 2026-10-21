"""M29 shared helpers: reuse the frozen common_ruler on the six original TRAIN contigs (CPU only).

Run with:  cd ~/coding-rna-reorg && PYTHONPATH=pipeline/m28:pipeline/m28/scripts:$HOME/envs/m29_extra \
           ~/envs/generanno/bin/python experiments/M29-CARBON-A-GAP/scripts/<script>.py
"""
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]          # ~/coding-rna-reorg
DATA = ROOT / 'data/m1_screen'
OUTS = ROOT / 'outputs'
M29 = OUTS / 'M29-CARBON-A-GAP'
CACHE = M29 / '_cache'
sys.path[:0] = [str(ROOT / 'pipeline/m28'), str(ROOT / 'pipeline/m28/scripts')]
import common_ruler as C   # noqa: E402  (frozen ruler: M26 CDS-assessable primary, exact CDS chain)

SCOPE = (('arabidopsis_thaliana', 'NC_003070.9', 30427671),
         ('arabidopsis_thaliana', 'NC_003071.7', 19698289),
         ('arabidopsis_thaliana', 'NC_003076.8', 26975502),
         ('oryza_sativa', 'NC_089043.1', 27474823),
         ('oryza_sativa', 'NC_089044.1', 23931887),
         ('oryza_sativa', 'NC_089046.1', 28271460))
DEV_SCOPE = (('arabidopsis_thaliana', 'NC_003074.8', 23459830),
             ('oryza_sativa', 'NC_089041.1', 29936421))
SEQ2SPECIES = {s: sp for sp, s, _ in SCOPE + DEV_SCOPE}
LENGTH = {s: n for _, s, n in SCOPE + DEV_SCOPE}
EXPECTED_PRIMARY = {'NC_003070.9': 7135, 'NC_003071.7': 4293, 'NC_003076.8': 6311,
                    'NC_089043.1': 1638, 'NC_089044.1': 1601, 'NC_089046.1': 1710}

# Per-tool TRAIN prediction GFFs (all six contigs), relative to OUTS.
TOOLS = {
    'CarbonA_unfiltered': 'M28-CARBON-A-TRAIN-RAW/{seq}/native/prediction.gff',
    'ANNEVO': 'M28-NATIVE-TRAIN-BASELINES/ANNEVO/{seq}/prediction.gff3',
    'Tiberius': 'M28-NATIVE-TRAIN-BASELINES/Tiberius/{seq}/prediction.gff3',
    'Helixer': 'M28-NATIVE-TRAIN-BASELINES/Helixer/{seq}/prediction.gff3',
    'B_freshjoint': 'M28-FULL-GRAMMAR-FRESH-JOINT24-TRAIN-RAW/{seq}/prediction/full_grammar.gff3',
    'E_lora': 'M28-FULL-GRAMMAR-JOINT-LORA-TRAIN-RAW/{seq}/prediction/full_grammar.gff3',
    'M_denseqv': 'M28-FULL-GRAMMAR-DENSE-QV-TRAIN-RAW/{seq}/prediction/full_grammar.gff3',
    'EM_pool': 'M28-FULL-GRAMMAR-E-M-ENERGY-POOL-TRAIN/{seq}/full_grammar.gff3',
}


def key(tx):
    return C.S.key(tx)           # (seqid, strand, ((a,b),...)); phase ignored


def gff_path(tool, seq):
    return OUTS / TOOLS[tool].format(seq=seq)


def read_sequence(seq):
    sp = SEQ2SPECIES[seq]
    return C.M.read_allowed(DATA / sp / 'genome.fa', {seq: LENGTH[seq]})[seq].upper()


def build_reference(seq):
    """Frozen ruler objects for one contig; cached as a pickle."""
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f'ref_{seq}.pkl'
    if p.exists():
        return pickle.loads(p.read_bytes())
    sp, length = SEQ2SPECIES[seq], LENGTH[seq]
    path = DATA / sp
    sequences = C.M.read_allowed(path / 'genome.fa', {seq: length})
    coding = C.S.E.parse_annotation(path / 'reference.gff3', {seq: length}, protein_coding_only=True)
    annotation = C.S.E.parse_annotation(path / 'reference.gff3', {seq: length})
    with (path / 'reference.gff3').open() as handle:
        _, _, raw_cds, _ = C.P.raw_features(handle, seq)
    refs, isoforms = C.select_reference(coding, raw_cds, sequences)
    if seq in EXPECTED_PRIMARY:
        assert len(refs) == EXPECTED_PRIMARY[seq], (seq, len(refs))
    genes = [(g['start'], g['end']) for g in annotation['genes'].values()]
    obj = dict(seq=seq, species=sp, length=length, refs=refs, isoforms=isoforms, gene_spans=genes,
               annotation=annotation)
    p.write_bytes(pickle.dumps(obj, protocol=4))
    return obj


def load_txs(path, seq, complete_only=False):
    pred = C.S.E.parse_annotation(Path(path), {seq: LENGTH[seq]})
    return C.S.E.primary_transcripts(pred, complete_only=complete_only)


def prf(tp, p, r):
    return C.S.prf(tp, p, r)
