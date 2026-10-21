"""M30 shared setup (protocol: docs/review_20261016/13_M30_protocol.md).  Imports the A1/A2 modules unchanged."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M29-A1-RESIDUAL' / 'scripts'))
sys.path.insert(0, str(HERE.parents[1] / 'M29-A2-TEMPERATURE' / 'scripts'))
import a1_common as K          # noqa: E402
OUT = K.OUTS / 'M30-HOLDOUT-ADAPT'
DATA = K.ROOT / 'data' / 'm1_screen'
RAW = K.ROOT / 'data' / 'm30' / 'raw'
VAL_BP = 20_000_000

SPECIES = {
    'dryas_octopetala': dict(
        asm='GCF_963921425.1_drDryOcto1.1',
        chroms=['NC_141366.1', 'NC_141367.1', 'NC_141368.1', 'NC_141369.1', 'NC_141370.1', 'NC_141371.1', 'NC_141372.1', 'NC_141373.1', 'NC_141374.1']),
    'eschscholzia_californica': dict(
        asm='GCF_041429995.1_dmEscCali1.0.p',
        chroms=['NW_028238198.1', 'NW_028238309.1', 'NW_028238420.1', 'NW_028238461.1', 'NW_028238472.1', 'NW_028238483.1']),
}
for sp in SPECIES.values():                     # protocol §2: sorted by accession, alternate A,B,A,B,...
    c = sorted(sp['chroms'])
    sp['A'] = c[0::2]; sp['B'] = c[1::2]
    sp['val_chrom'] = sp['A'][-1]
    sp['val_contig'] = sp['val_chrom'] + '__val20M'
    sp['train'] = sp['A'][:-1]


def contigs(species):
    s = SPECIES[species]
    return dict(train=list(s['train']), val=[s['val_contig']], B=list(s['B']))


def lengths(species):
    p = DATA / species / 'lengths.json'
    if p.exists():
        return {k: int(v) for k, v in json.load(open(p)).items()}
    import m29_common as M29          # legacy species (arabidopsis/rice): lengths known to the M29 helpers
    return {s: M29.LENGTH[s] for s, sp in M29.SEQ2SPECIES.items() if sp == species}


def register():
    """Make the M29 helpers (reference builder, sequence reader, GFF loader) aware of the M30 contigs."""
    import m29_common as M29
    for sp in SPECIES:
        p = DATA / sp / 'lengths.json'
        if p.exists():
            for k, v in json.load(open(p)).items():
                M29.SEQ2SPECIES[k] = sp
                M29.LENGTH[k] = int(v)
    return M29


def species_of(contig):
    for sp in SPECIES:
        if contig in lengths(sp):
            return sp
    raise KeyError(contig)


def raw_parquet(contig):
    import glob
    sp = species_of(contig)
    p = sorted(glob.glob(str(OUT / 'carbonA_raw' / sp / contig / '*.raw_probs.parquet')))
    assert len(p) == 1, (contig, p)
    return p[0]


def offset_path(contig, strand):
    return OUT / 'offset' / f'{contig}_{strand}.npy'
