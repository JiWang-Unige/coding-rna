"""M34-NTV2-SCALEUP shared setup (protocol: docs/review_20261016/19_M34_scaleup_plan.md).
Imports the frozen M28/M29 code unchanged (never edits it).  CPU-only data preparation; no GPU, no model.

Run on ont:  cd ~/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts && taskset -c 24-27 nice -n 10 ~/envs/generanno/bin/python <script>.py
"""
import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if 'a1_common' not in sys.modules:       # (the M32 PlantCAD2 shim may already be loaded: py3.11 env, must stay first on the path)
    sys.path.insert(0, str(HERE.parents[1] / 'M29-A1-RESIDUAL' / 'scripts'))
import a1_common as K          # noqa: E402  (sets sys.path for pipeline/m28 and the M29 helpers)

ROOT = K.ROOT                                             # ~/coding-rna-reorg
DATA = ROOT / 'data' / 'm34_train'                        # -> ~/coding-rna/data/m34_train
RAW = DATA / 'raw'                                        # NCBI downloads (md5 verified)
PREP = DATA / 'prepared'                                  # per new species: genome.fa, reference.gff3, lengths.json, splits.json
CATALOG = DATA / 'catalog'                                # per species (all 10): transcripts.jsonl, windows_*.jsonl, summary.json, labeled_k*.json
OUT = K.OUTS / 'M34-NTV2-SCALEUP'
M1 = ROOT / 'data' / 'm1_screen'

VAL_BP = 8_000_000            # internal-validation crop of the held-out chromosome of each NEW species
MIN_CHROM_BP = 5_000_000      # assembled-molecule chromosomes shorter than this are dropped (Dioscorea has 8 sub-5-Mb pseudo-chromosomes)
SELECTION_SEED = 20261017     # same nested labeled-gene permutation seed as the M30 protocol (13 §3)
WIDTH, STRIDE, FLANK = 24576, 12288, 256

# ---- the two legacy species (already used by M28/M31): training chromosomes = the four M31 TRAIN chromosomes -----------------------
LEGACY = {
    'arabidopsis_thaliana': dict(order=0, organism='Arabidopsis thaliana', accession='GCF_000001735.4', data=M1 / 'arabidopsis_thaliana',
                                 train=['NC_003070.9', 'NC_003076.8'], val=['NC_003071.7'], val_crop=None,
                                 never_train=['NC_003074.8', 'NC_003075.7'], family='Brassicaceae'),
    'oryza_sativa': dict(order=0, organism='Oryza sativa', accession='GCF_034140825.1', data=M1 / 'oryza_sativa',
                         train=['NC_089043.1', 'NC_089046.1'], val=['NC_089044.1'], val_crop=None,
                         never_train=['NC_089041.1', 'NC_089042.1'], family='Poaceae'),
}
# DEV chromosomes (evaluated ONCE per final arm) and never-used chromosomes of the legacy species; the sampler asserts none is ever sampled.
DEV_CONTIGS = ('NC_003074.8', 'NC_089041.1')

FAMILY = {'populus_trichocarpa': 'Salicaceae', 'brachypodium_distachyon': 'Poaceae', 'solanum_lycopersicum': 'Solanaceae',
          'phoenix_dactylifera': 'Arecaceae', 'phaseolus_vulgaris': 'Fabaceae', 'citrus_sinensis': 'Rutaceae',
          'daucus_carota': 'Apiaceae', 'dioscorea_rotundata': 'Dioscoreaceae'}


def new_species():
    rows = json.load(open(HERE / 'species.json'))
    return {r['key']: r for r in sorted(rows, key=lambda r: r['order'])}


def raw_dir(row):
    return RAW / f"{row['accession']}_{row['assembly_name']}"


def parse_assembly_report(path):
    """-> list of dicts for every sequence of the NCBI assembly report."""
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        rows.append(dict(name=f[0], role=f[1], molecule=f[2], location=f[3], genbank=f[4], refseq=f[6], unit=f[7], length=int(f[8])))
    return rows


def main_chromosomes(report_rows):
    """Project rule (main nuclear chromosomes only): role assembled-molecule, location type Chromosome, Primary Assembly, length >= MIN_CHROM_BP."""
    return [r for r in report_rows if r['role'] == 'assembled-molecule' and r['location'] == 'Chromosome'
            and r['unit'] == 'Primary Assembly' and r['length'] >= MIN_CHROM_BP]


def choose_heldout(chroms):
    """Deterministic, length-only rule: the lower-median-length chromosome (ties by RefSeq accession) is the internal-validation chromosome
    of a NEW species; only its first VAL_BP bases are used (genes fully inside), the rest is discarded (neither training nor evaluation)."""
    order = sorted(chroms, key=lambda r: (r['length'], r['refseq']))
    return order[(len(order) - 1) // 2]


def species_table():
    """key -> dict(train=[contigs], val=[contig names], lengths={contig: len}, dir=Path to genome.fa/reference.gff3)."""
    table = {}
    for key, sp in LEGACY.items():
        lengths = {}
        with open(sp['data'] / 'genome.fa') as f:                       # headers + sequence lengths of the wanted contigs only
            q, n = None, 0
            for line in f:
                if line[0] == '>':
                    if q in sp['train'] + sp['val']:
                        lengths[q] = n
                    q, n = line[1:].split()[0], 0
                else:
                    n += len(line.strip())
            if q in sp['train'] + sp['val']:
                lengths[q] = n
        table[key] = dict(sp, key=key, dir=sp['data'], lengths=lengths)
    for key, row in new_species().items():
        p = PREP / key / 'splits.json'
        if p.exists():
            s = json.load(open(p))
            table[key] = dict(key=key, organism=row['organism'], accession=row['accession'], order=row['order'], dir=PREP / key,
                              train=s['train'], val=[s['validation_contig']], lengths={k: int(v) for k, v in s['lengths'].items()},
                              val_crop=s['validation_crop'], family=FAMILY[key])
    return table


def select_labeled(rows, k, seed=SELECTION_SEED):
    """Nested labeled-gene selection (M30 protocol 13 §3): rows = eligible transcripts (dicts with seqid, span, id); first ceil(k*n) of a fixed permutation of the
    list sorted by (seqid, span start, id).  k1 < k2  =>  labeled(k1) is a subset of labeled(k2)."""
    import math
    import numpy as np
    elig = sorted(rows, key=lambda t: (t['seqid'], t['span'][0], t['id']))
    perm = np.random.default_rng(seed).permutation(len(elig))
    n = int(math.ceil(k * len(elig) - 1e-9))
    return {elig[i]['id'] for i in perm[:n]}


_GRID = {}


def grid(length):
    if length not in _GRID:
        import manifest as MF
        cells = [(a, min(length, a + WIDTH)) for a in MF.starts(length)]
        _GRID[length] = (cells, [a for a, _ in cells])
    return _GRID[length][0]


def gene_cells(span, length):
    """Grid windows (a, b) overlapping the CDS span [s, e) — the same cells as `Labeler.labeled_windows` of M30 (±2 neighbours around the left-most)."""
    import bisect
    cells = grid(length)
    starts = _GRID[length][1]
    i = bisect.bisect_right(starts, span[0]) - 1
    return [(a, b) for a, b in cells[max(0, i - 2):min(len(cells), i + 3)] if a < span[1] and b > span[0]]
