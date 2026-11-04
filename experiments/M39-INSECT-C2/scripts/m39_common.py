"""M39-INSECT-C2 shared setup (protocol: docs/review_20261016/29_animal_C2_insect_protocol.md).  Copy of the M34 `m34_common.py` with the insect species:
every species (2 fully labeled old species + 8 new species at 50 %) is prepared from the raw NCBI files into data/m39_insect/prepared/<key>/.
Imports the frozen M28/M29 code unchanged (never edits it).  CPU-only data preparation; no GPU, no model.

Run on ont:  cd ~/coding-rna-reorg/experiments/M39-INSECT-C2/scripts && CUDA_VISIBLE_DEVICES= taskset -c 24-31 nice -n 19 ~/envs/generanno/bin/python <script>.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if 'a1_common' not in sys.modules:
    sys.path.insert(0, str(HERE.parents[1] / 'M29-A1-RESIDUAL' / 'scripts'))
import a1_common as K          # noqa: E402  (sets sys.path for pipeline/m28 and the M29 helpers)

ROOT = K.ROOT                                             # ~/coding-rna-reorg
DATA = ROOT / 'data' / 'm39_insect'                       # -> ~/coding-rna/data/m39_insect
RAW = DATA / 'raw'                                        # NCBI downloads (md5 verified)
PREP = DATA / 'prepared'                                  # per species: genome.fa, reference.gff3, lengths.json, splits.json
CATALOG = DATA / 'catalog'                                # per species: transcripts.jsonl, windows_*.jsonl, summary.json, labeled_k*.json
OUT = K.OUTS / 'M39-INSECT-C2'

MIN_CHROM_BP = 5_000_000      # assembled-molecule chromosomes shorter than this are dropped (as M34)
LOCATIONS = ('Chromosome', 'Linkage Group')               # protocol 29 §2.1 E1: Apis cerana / Ips typographus register their chromosomes as linkage groups
VAL_MIN_BP = 8_000_000        # validation crop of a NEW species: at least 8 Mb ...
VAL_STEP_BP = 1_000_000       # ... grown by whole Mb ...
VAL_TARGET_REFS = 500         # ... until it holds >= 500 CDS-assessable primary references (genes wholly inside), capped at the whole chromosome (protocol 29 §3.2)
VAL_MIN_REFS = 400            # E8: fewer than this in the final crop -> replace the species (protocol 29 §2.4)
SELECTION_SEED = 20261017     # same nested labeled-gene permutation seed as M34 / the M30 protocol (13 §3)
WIDTH, STRIDE, FLANK = 24576, 12288, 256


def species_rows():
    return json.load(open(HERE / 'species.json'))


def roster():
    """key -> species.json row for the species in the experiment (role old / new), in protocol order."""
    return {r['key']: r for r in sorted(species_rows(), key=lambda r: r['order']) if r['role'] in ('old', 'new')}


def old_keys():
    return [k for k, r in roster().items() if r['role'] == 'old']


def new_keys():
    return [k for k, r in roster().items() if r['role'] == 'new']


def raw_dir(row):
    return RAW / row['ftp_stem']


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
    """Project rule (main nuclear chromosomes only): role assembled-molecule, location Chromosome or Linkage Group, Primary Assembly, length >= MIN_CHROM_BP,
    with a RefSeq accession (the Apis cerana report also lists the GenBank-only copies of its linkage groups, RefSeq-Accn 'na')."""
    return [r for r in report_rows if r['role'] == 'assembled-molecule' and r['location'] in LOCATIONS
            and r['unit'] == 'Primary Assembly' and r['length'] >= MIN_CHROM_BP and r['refseq'] != 'na']


def length_order(chroms):
    return sorted(chroms, key=lambda r: (r['length'], r['refseq']))


def choose_heldout(chroms):
    """Lower-median-length chromosome (ties by RefSeq accession): the validation chromosome of a NEW species / the selection chromosome S of an OLD species."""
    order = length_order(chroms)
    return order[(len(order) - 1) // 2]


def choose_d(chroms):
    """OLD species only: the descriptive held-out chromosome D = the next shorter one in the same length order (protocol 29 §3.2)."""
    order = length_order(chroms)
    return order[(len(order) - 1) // 2 - 1]


def species_table():
    """key -> dict(role, train=[contigs], val=[selection chromosome S (old) | validation crop (new)], heldout_d=[D (old) | none], lengths, dir, family)."""
    table = {}
    for key, row in roster().items():
        p = PREP / key / 'splits.json'
        if p.exists():
            s = json.load(open(p))
            table[key] = dict(key=key, role=row['role'], organism=row['organism'], accession=row['accession'], order=row['order'], dir=PREP / key,
                              train=s['train'], val=[s['validation_contig']], heldout_d=list(s.get('heldout_d', [])),
                              lengths={k: int(v) for k, v in s['lengths'].items()}, val_crop=s.get('validation_crop'), family=row['family'])
    return table


def forbidden_contigs():
    """Contigs that must never be sampled: every S / validation crop and every D."""
    out = set()
    for v in species_table().values():
        out |= set(v['val']) | set(v['heldout_d'])
    return out


def select_labeled(rows, k, seed=SELECTION_SEED):
    """Nested labeled-gene selection (M30 protocol 13 §3, same as M34): rows = eligible transcripts (dicts with seqid, span, id); first ceil(k*n) of a fixed
    permutation of the list sorted by (seqid, span start, id).  k1 < k2  =>  labeled(k1) is a subset of labeled(k2)."""
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
    """Grid windows (a, b) overlapping the CDS span [s, e) -- the same cells as `Labeler.labeled_windows` of M30 (+-2 neighbours around the left-most)."""
    import bisect
    cells = grid(length)
    starts = _GRID[length][1]
    i = bisect.bisect_right(starts, span[0]) - 1
    return [(a, b) for a, b in cells[max(0, i - 2):min(len(cells), i + 3)] if a < span[1] and b > span[0]]
