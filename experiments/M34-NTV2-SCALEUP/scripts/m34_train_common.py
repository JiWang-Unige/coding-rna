"""M34 shared training/inference helpers: multi-species labels (from the catalogs written by m34_labels.py, no GFF re-parse), the draw -> example
builder, and the registration of the NEW-species validation pseudo-contigs in the frozen M29 helpers (reference builder / sequence reader), so the
unchanged M31 inference (`c1_infer`) and validation scoring (`c1_val`) code paths work on them.  Frozen code and M31/M32/M33 files are only imported.
"""
import json
import pickle
import re
import sys
from pathlib import Path

import os                                   # noqa: E402

HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
BACKEND = os.environ.get('M34_BACKEND', 'ntv2')          # 'ntv2' (NT-v2-500M + q/v LoRA, M31 B recipe; env generanno) | 'pc2' (PlantCAD2-Small, M32 geometry; env plantcad2)
assert BACKEND in ('ntv2', 'pc2'), BACKEND
ENC_SCRIPTS = EXPS / ('M31-C1-NTV2' if BACKEND == 'ntv2' else 'M32-C1-PLANTCAD2') / 'scripts'
for p in (str(ENC_SCRIPTS), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np                          # noqa: E402

if BACKEND == 'ntv2':
    import ntv2_adapter as A                # noqa: E402,F401  (M31; sets the pipeline/M29 sys.path)
    import m31_common as M                  # noqa: E402  (G, Q, lr_at, atomic_save, parse_cpus)
else:
    import pc2_adapter as A                 # noqa: E402,F401  (M32; sets the pipeline/M29 sys.path, M32 `a1_common` shim first -- py3.11 plantcad2 env)
    import m32_common as M                  # noqa: E402  (same helpers as m31_common)
import c1_model as MD                       # noqa: E402  (the backend's head/feature plumbing: M31 or M32 version, resolved by the path above)
import m34_common as C                      # noqa: E402
import m34_labels as L                      # noqa: E402
import m34_schedule as S                    # noqa: E402

G, Q = M.G, M.Q
MF = L.MF
if BACKEND == 'ntv2':
    INIT = M.OUT / 'runA' / 'model_002048.pt'           # M31 training-A head (frozen encoder, step 2048): the common initialisation of every arm
else:
    _best = M.OUT / 'runA' / 'best_step.txt'            # M32 frozen-PlantCAD2 head selected by its internal validation
    INIT = M.OUT / 'runA' / f"model_{int(_best.read_text()):06d}.pt" if _best.exists() else M.OUT / 'runA' / 'model_002048.pt'
OUT = C.OUT
LEGACY = tuple(S.LEGACY)
VAL_REGIONS_LEGACY = ('NC_003071.7', 'NC_089044.1')     # whole held-out chromosomes (selection set of every arm)
SCHEDULES = OUT / 'schedules'
CACHE = OUT / '_cache'


def arm_dir(arm):
    """Output directory of an arm.  NT-v2 keeps the protocol path <OUT>/<arm>; the PlantCAD2 variant lives next to it as <OUT>/pc2_<arm>."""
    return OUT / (arm if BACKEND == 'ntv2' else f'pc2_{arm}')


# ------------------------------------------------------------------------------------------------------------------------------
def new_keys():
    return [k for k in C.species_table() if k not in LEGACY]


def val_contig(key):
    return C.species_table()[key]['val'][0]


def arm_spec(arm):
    return S.ARMS[arm]


def arm_val_regions(arm):
    """-> (selection regions [Ath/Osa held-out chromosomes], arm's own new-species validation pseudo-contigs)"""
    own = [val_contig(k) for k in S.ARMS[arm]['species'] if k not in LEGACY]
    return list(VAL_REGIONS_LEGACY), own


def other_val_regions(arm):
    """new-species validation pseudo-contigs of the species NOT in the arm (cross-species transfer set)"""
    have = set(S.ARMS[arm]['species'])
    return [val_contig(k) for k in new_keys() if k not in have]


# ---- registration of the pseudo-contigs in the frozen M29 helpers (no frozen file is edited) -------------------------------
_REG = {}


def register():
    """Make m29_common (build_reference / read_sequence / load_txs / LENGTH / SEQ2SPECIES) aware of the M34 validation pseudo-contigs.
    Idempotent; must also run in every spawned decode worker (see m34_infer._decode_init)."""
    import m29_common as MC
    if _REG:
        return MC
    tab = C.species_table()
    for key in new_keys():
        sp = tab[key]
        for contig in sp['val']:
            MC.SEQ2SPECIES[contig] = key
            MC.LENGTH[contig] = int(sp['lengths'][contig])
            _REG[contig] = key
    orig_read, orig_build = MC.read_sequence, MC.build_reference
    seq_cache = {}

    def read_sequence(seq):
        if seq not in _REG:
            return orig_read(seq)
        if seq not in seq_cache:
            seq_cache.clear()
            seq_cache[seq] = MC.C.M.read_allowed(C.PREP / _REG[seq] / 'genome.fa', {seq: MC.LENGTH[seq]})[seq].upper()
        return seq_cache[seq]

    def build_reference(seq):
        """Same construction as the frozen `m29_common.build_reference` for a pseudo-contig of a prepared (m34) species."""
        if seq not in _REG:
            return orig_build(seq)
        CACHE.mkdir(parents=True, exist_ok=True)
        pk = CACHE / f'ref_{seq}.pkl'
        if pk.exists():
            return pickle.loads(pk.read_bytes())
        key, length = _REG[seq], MC.LENGTH[seq]
        path = C.PREP / key
        sequences = MC.C.M.read_allowed(path / 'genome.fa', {seq: length})
        coding = MC.C.S.E.parse_annotation(path / 'reference.gff3', {seq: length}, protein_coding_only=True)
        annotation = MC.C.S.E.parse_annotation(path / 'reference.gff3', {seq: length})
        with (path / 'reference.gff3').open() as handle:
            _, _, raw_cds, _ = MC.C.P.raw_features(handle, seq)
        refs, isoforms = MC.C.select_reference(coding, raw_cds, sequences)
        genes = [(g['start'], g['end']) for g in annotation['genes'].values()]
        obj = dict(seq=seq, species=key, length=length, refs=refs, isoforms=isoforms, gene_spans=genes, annotation=annotation)
        pk.write_bytes(pickle.dumps(obj, protocol=4))
        return obj

    MC.read_sequence, MC.build_reference = read_sequence, build_reference
    return MC


# ---- labels ---------------------------------------------------------------------------------------------------------------
class CatalogLabeler(L.Labeler):
    """`m34_labels.Labeler` restored from the catalog (`transcripts.jsonl`) instead of re-parsing the GFF3 (minutes per species).  Same `targets`, `window`,
    `local` code as the parent; the catalog rows are exactly what the parent's constructor built.  Only training contigs are loaded."""

    def __init__(self, key, contigs, k=1.0, selection_seed=C.SELECTION_SEED):
        sp = C.species_table()[key]
        self.species, self.contigs, self.k = key, list(contigs), k
        self.lengths = {c: sp['lengths'][c] for c in self.contigs}
        self.seqs = MF.read_allowed(sp['dir'] / 'genome.fa', self.lengths)
        keep = set(self.contigs)
        self.catalog = {}
        with open(C.CATALOG / key / 'transcripts.jsonl') as f:
            for line in f:
                t = json.loads(line)
                if t['seqid'] in keep:
                    t['CDS'] = [tuple(x) for x in t['CDS']]
                    self.catalog[t['id']] = t
        elig = [t for t in self.catalog.values() if t['eligible']]
        self.labeled = C.select_labeled(elig, k, selection_seed)
        self.n_eligible, self.n_labeled = len(elig), len(self.labeled)
        from collections import defaultdict
        self.by_qs = defaultdict(list)
        for t in self.catalog.values():
            self.by_qs[t['seqid'], t['strand']].append(t)
        self.starts, self.max_span = {}, {}
        for key2, lst in self.by_qs.items():
            lst.sort(key=lambda t: t['span'][0])
            self.starts[key2] = [t['span'][0] for t in lst]
            self.max_span[key2] = max((t['span'][1] - t['span'][0] for t in lst), default=0)


def arm_labelers(arm, species=None):
    """{species: CatalogLabeler over its TRAINING contigs at the arm's label fraction}, treatment name.  Legacy species are always fully labeled."""
    spec = S.ARMS[arm]
    tab = C.species_table()
    labs = {}
    for key in (species or spec['species']):
        k = 1.0 if key in LEGACY else spec['new_k']
        labs[key] = CatalogLabeler(key, tab[key]['train'], k)
    return labs, spec['treatment']


def example(lab, treatment, row):
    """Draw row -> training example of the frozen geometry (24,576-bp core +- one 12,288 neighbour, oriented context, core targets).
    Same construction as the M31 `m31_common.example` / M30 `example`, but the targets come from the M34 labeler (FULL / U / F)."""
    c, s, a, b = row['seqid'], row['strand'], row['start'], row['end']
    tg, win, core = lab.targets(c, s, a, b, treatment)
    chromosome = lab.seqs[c]
    geometry = G.context_geometry(a, b, len(chromosome))
    ca, cb = geometry['genomic_context']
    n = b - a
    seg, region = tg['seg'][:n], tg['region'][:n]
    fine, coarse = int((seg >= 0).sum()), int(((seg < 0) & (region >= 0)).sum())
    offset = a - ca if s == '+' else cb - b
    context = re.sub('[^ACGT]', 'N', chromosome[ca:cb])
    if s == '-':
        context = MF.rc(context)
    return chromosome, geometry, context, offset, dict(seg=seg, region=region), (fine, coarse, n - fine - coarse), win


def read_schedule(path):
    with open(path) as f:
        return [json.loads(x) for x in f]


def forbidden_contigs():
    """Contigs that must never be sampled: DEV, the never-used legacy chromosomes, every validation chromosome / pseudo-contig."""
    tab = C.species_table()
    out = set(C.DEV_CONTIGS)
    for s in C.LEGACY.values():
        out |= set(s['val']) | set(s['never_train'])
    for v in tab.values():
        out |= set(v['val'])
    return out
