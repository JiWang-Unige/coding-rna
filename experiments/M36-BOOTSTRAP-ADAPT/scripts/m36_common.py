"""M36-BOOTSTRAP-ADAPT shared setup (protocol: docs/review_20261016/22_M36_bootstrap_adaptation_protocol.md, CPU preparation §7).
Imports the M30 / M34 helpers and the frozen M28 code unchanged.  CPU only: no model, no GPU, no target-species reference annotation.

Run on ont:  cd ~/coding-rna-reorg/experiments/M36-BOOTSTRAP-ADAPT/scripts && taskset -c 24-31 nice -n 10 ~/envs/generanno/bin/python <script>.py ...
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
for p in (str(EXPS / 'M34-NTV2-SCALEUP' / 'scripts'), str(EXPS / 'M30-HOLDOUT-ADAPT' / 'scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)
import m34_common as C          # noqa: E402  (sets the M29-A1 / pipeline/m28 sys.path through a1_common)
import m30_common as M30        # noqa: E402  (SPECIES, A/B split, DATA = data/m1_screen)
import a1_common as K           # noqa: E402

ROOT = C.ROOT                                            # ~/coding-rna-reorg
DATA = ROOT / 'data' / 'm36'                             # -> ~/coding-rna/data/m36   (proteins, A-only genomes, miniprot raw, seeds)
OUT = K.OUTS / 'M36-BOOTSTRAP-ADAPT'                     # -> ~/coding-rna/outputs/M36-BOOTSTRAP-ADAPT
M1 = M30.DATA                                            # data/m1_screen/<species>/{genome.fa, inputs/<contig>.fasta, lengths.json}
MINIPROT = os.path.expanduser('~/envs/miniprot/bin/miniprot')

# ---- S2 seed filter (protocol 22 §2, fixed before any alignment was run) --------------------------------------------------------------
COVERAGE_MIN = 0.95        # aligned query span / protein length
IDENTITY_MIN = 0.70        # miniprot Identity of the alignment
MAX_EXONS = 128            # frozen grammar exon cap (the labels/DP use 128)
PSEUDO_BACKGROUND_MIN_P = 0.99   # S1/S3 pseudo-background: per-base P(no gene covers the base) >= this

# ---- target species: label-set A contigs only (M30 protocol 13 §2).  Evaluation-set B contigs are NEVER listed here. --------------------
TARGETS = {}
for _sp, _row in M30.SPECIES.items():
    TARGETS[_sp] = dict(asm=_row['asm'], train=list(_row['train']), val=[_row['val_contig']], A=list(_row['train']) + [_row['val_contig']],
                        B_never_touched=list(_row['B']))


def a_lengths(species):
    """{A contig: length} from data/m1_screen/<species>/lengths.json (train chromosomes + the 20-Mb validation pseudo-contig)."""
    lens = M30.lengths(species)
    return {c: lens[c] for c in TARGETS[species]['A']}


def assert_not_B(species, contigs):
    bad = set(contigs) & set(TARGETS[species]['B_never_touched'])
    if bad:
        raise RuntimeError(f'evaluation-set B contig(s) requested for {species}: {sorted(bad)}')


# ---- protein sources: the 10 M34 training species (RefSeq proteins already on ont, data/m34_train/raw) ----------------------------------
def protein_sources():
    """-> [dict(key, organism, faa, gff, keep_contigs|None)] in the fixed M34 species order (legacy first)."""
    rows = []
    for key, sp in C.LEGACY.items():
        acc = {'arabidopsis_thaliana': 'GCF_000001735.4_TAIR10.1', 'oryza_sativa': 'GCF_034140825.1_ASM3414082v1'}[key]
        rows.append(dict(key=key, organism=sp['organism'], faa=C.RAW / 'legacy_proteins' / f'{acc}_protein.faa.gz',
                         gff=sp['data'] / 'reference.gff3', keep_contigs=set(sp['train']) | set(sp['val'])))   # never the DEV / split_test chromosomes
    for key, row in C.new_species().items():
        d = C.raw_dir(row)
        name = d.name
        rows.append(dict(key=key, organism=row['organism'], faa=d / f'{name}_protein.faa.gz', gff=d / f'{name}_genomic.gff.gz', keep_contigs=None))
    return rows


def pin_and_nice(cpus=range(24, 32), niceness=10):
    try:
        os.sched_setaffinity(0, set(cpus))
    except (AttributeError, OSError):
        pass
    try:
        os.nice(niceness)
    except OSError:
        pass
