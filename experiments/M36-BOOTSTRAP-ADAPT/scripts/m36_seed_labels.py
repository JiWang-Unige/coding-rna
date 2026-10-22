"""Seed GFF -> M34 label interface for the M36 arms (protocol 22 §2/§3).  Drop-in for `m34_labels.Labeler` in `m34_train_common.example`
(same `targets(contig, strand, a, b, treatment) -> (dict(seg, region), win, core)`, `seqs`, `lengths`, `local`, `window`, `grid`).
  S1 / S2 / S3 : seed loci get the fine labels of the frozen `build_targets` (seed = complete supported primary transcript); everything else is
                 unknown (IGNORE in both the fine `seg` and the coarse `region` channel, i.e. marginalised), except optional pseudo-background
                 intervals (S1/S3: per-base P(no gene) >= 0.99 segments from m36_nogene_marginal) which are labelled background (0 / 0).
  F1           : seeds labelled, every other position background (the Carbon-A-style fill-background control; no IGNORE except non-ACGT / padding).
Seeds come from a seed GFF3 (m36_seed_filter output, or any GFF3 with mRNA/CDS lines whose CDS includes the stop codon) or from seed rows (JSON).
No reference annotation of the target species is read here.
"""
import bisect
import json
from collections import defaultdict

import numpy as np

import m36_common as W
import m34_common as C
import manifest as MF                                   # frozen: read_allowed, orient_interval, rc, merged
from src.m28.labels import build_targets, IGNORE        # frozen

WIDTH = MF.WIDTH
S_ARMS = ('S', 'S1', 'S2', 'S3')
F_ARMS = ('F', 'F1', 'FULL')


# ---- seed rows ----------------------------------------------------------------------------------------------------------------------
def seed_row(seqid, strand, cds, ident, **extra):
    """cds: [(a0, b0, phase)] ascending, stop codon included; -> catalog row with the keys `build_targets` / `Labeler.targets` use."""
    cds = [(int(a), int(b), str(p)) for a, b, p in sorted(cds)]
    row = dict(id=ident, gene_id=extra.pop('gene_id', ident), seqid=seqid, strand=strand, CDS=cds, span=[cds[0][0], cds[-1][1]],
               training_primary=True, CDS_complete_supported=True, overlapping_primary=False, assessable=True, eligible=True, eligible_positive=True,
               n_exons=len(cds), cds_bp=sum(b - a for a, b, _ in cds))
    row.update(extra)
    return row


def phases(cds_asc, strand):
    order = cds_asc if strand == '+' else cds_asc[::-1]
    out, off = {}, 0
    for a, b in order:
        out[(a, b)] = str((3 - off % 3) % 3)
        off += b - a
    return [(a, b, out[(a, b)]) for a, b in cds_asc]


def load_seed_gff(path, contigs=None):
    """GFF3 with mRNA (ID, Parent) + CDS (Parent) lines -> [seed rows]; coordinates 1-based inclusive -> 0-based half-open; phases recomputed."""
    def attrs(s):
        return dict(x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)
    mrna, cds = {}, defaultdict(list)
    with open(path) as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or (contigs is not None and p[0] not in contigs):
                continue
            a = attrs(p[8])
            if p[2] == 'mRNA':
                mrna[a['ID']] = dict(seqid=p[0], strand=p[6], parent=a.get('Parent', a['ID']), attrs=a)
            elif p[2] == 'CDS':
                cds[a['Parent']].append((int(p[3]) - 1, int(p[4])))
    rows = []
    for mid, m in mrna.items():
        if not cds.get(mid):
            continue
        ex = sorted(cds[mid])
        a = m['attrs']
        rows.append(seed_row(m['seqid'], m['strand'], phases(ex, m['strand']), a.get('Parent', mid), gene_id=m['parent'],
                             source_protein=a.get('source_protein'), source_species=a.get('source_species'),
                             identity=float(a['identity']) if 'identity' in a else None, coverage=float(a['coverage']) if 'coverage' in a else None))
    rows.sort(key=lambda r: (r['seqid'], r['span'][0], r['strand']))
    return rows


def load_seed_json(path, contigs=None):
    rows = json.load(open(path))
    return [seed_row(r['seqid'], r['strand'], r['CDS'], r['id'], **{k: v for k, v in r.items() if k not in ('seqid', 'strand', 'CDS', 'id', 'span', 'n_exons', 'cds_bp')})
            for r in rows if contigs is None or r['seqid'] in contigs]


def merge_seed_sets(primary, secondary):
    """S3 = S1 ∪ S2 with S2 priority: secondary seeds overlapping (same contig & strand) a primary seed are dropped.  -> (rows, n_dropped)"""
    idx = defaultdict(list)
    for r in primary:
        idx[r['seqid'], r['strand']].append(tuple(r['span']))
    for v in idx.values():
        v.sort()
    out, dropped = list(primary), 0
    for r in secondary:
        spans = idx.get((r['seqid'], r['strand']), [])
        a, b = r['span']
        i = bisect.bisect_left(spans, (b, b))
        if any(x < b and y > a for x, y in spans[max(0, i - 64):i + 1]):
            dropped += 1
            continue
        out.append(r)
    out.sort(key=lambda r: (r['seqid'], r['span'][0], r['strand']))
    return out, dropped


def pseudo_background(p_nogene_oriented, strand, length, threshold=W.PSEUDO_BACKGROUND_MIN_P, min_len=1):
    """per-base P(no gene covers the base) on one oriented strand -> genomic 0-based half-open intervals with P >= threshold (runs >= min_len)."""
    p = np.asarray(p_nogene_oriented)
    assert p.shape == (length,), (p.shape, length)
    hit = p >= threshold
    if not hit.any():
        return []
    d = np.diff(np.concatenate([[0], hit.astype(np.int8), [0]]))
    starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    iv = [(int(x), int(y)) for x, y in zip(starts, ends) if y - x >= min_len]
    if strand == '-':
        iv = [(length - y, length - x) for x, y in reversed(iv)]
    return iv


# ---- labeler ------------------------------------------------------------------------------------------------------------------------
class SeedLabeler:
    def __init__(self, species, genome_fa, lengths, seeds, background=None):
        """lengths: {contig: length} (label-set A contigs only); seeds: [seed rows] or a path (.gff3 / .json);
        background: {(contig, strand): [(a, b)]} genomic 0-based half-open pseudo-background intervals (S1/S3), optional."""
        self.species, self.contigs = species, list(lengths)
        self.lengths = dict(lengths)
        if species in W.TARGETS:
            W.assert_not_B(species, self.contigs)
        self.seqs = MF.read_allowed(genome_fa, self.lengths) if genome_fa is not None else None
        if isinstance(seeds, (str, bytes)) or hasattr(seeds, 'suffix'):
            seeds = load_seed_gff(seeds, set(self.contigs)) if str(seeds).endswith(('.gff3', '.gff')) else load_seed_json(seeds, set(self.contigs))
        self.catalog = {r['id']: r for r in seeds if r['seqid'] in self.lengths}
        self.labeled = set(self.catalog)
        self.n_eligible = self.n_labeled = len(self.catalog)
        self.by_qs = defaultdict(list)
        for t in self.catalog.values():
            self.by_qs[t['seqid'], t['strand']].append(t)
        self.starts, self.max_span = {}, {}
        for key, lst in self.by_qs.items():
            lst.sort(key=lambda t: t['span'][0])
            self.starts[key] = [t['span'][0] for t in lst]
            self.max_span[key] = max((t['span'][1] - t['span'][0] for t in lst), default=0)
        self.background = {}
        for key, iv in (background or {}).items():
            m = MF.merged([list(x) for x in iv])
            self.background[key] = ([a for a, _ in m], [b for _, b in m])

    # -- same geometry helpers as m34_labels.Labeler --
    def grid(self, contig):
        return C.grid(self.lengths[contig])

    def window(self, contig, strand, a, b):
        return dict(id=f'{self.species}|{contig}|{a}|{strand}', species=self.species, seqid=contig, start=a, end=b, width=WIDTH,
                    valid_bases=b - a, strand=strand)

    def local(self, contig, strand, a, b):
        key = (contig, strand)
        lst = self.by_qs.get(key, [])
        if not lst:
            return []
        lo = bisect.bisect_left(self.starts[key], a - self.max_span[key])
        hi = bisect.bisect_left(self.starts[key], b)
        return [t for t in lst[lo:hi] if t['span'][1] > a and t['span'][0] < b]

    def background_local(self, contig, strand, a, b):
        key = (contig, strand) if (contig, strand) in self.background else contig
        if key not in self.background:
            return []
        s, e = self.background[key]
        i = bisect.bisect_left(e, a + 1)
        out = []
        while i < len(s) and s[i] < b:
            out.append((max(a, s[i]), min(b, e[i])))
            i += 1
        return out

    def seed_cells(self):
        """{seed id: [(contig, strand, a, b)]} grid cells covering each seed (the 'seed-first' half of the window sampler)."""
        return {tid: [(t['seqid'], t['strand'], a, b) for a, b in C.gene_cells(t['span'], self.lengths[t['seqid']])] for tid, t in self.catalog.items()}

    def targets(self, contig, strand, a, b, treatment):
        """treatment: 'S1' | 'S2' | 'S3' | 'S' (seed loci labelled, rest unknown, optional pseudo-background) | 'F1' | 'F' | 'FULL' (rest background)."""
        if treatment not in S_ARMS + F_ARMS:
            raise ValueError(treatment)
        local = self.local(contig, strand, a, b)
        win = self.window(contig, strand, a, b)
        win.update(overlapping_transcript_ids=[t['id'] for t in local], positive_ids=[], chain_unknown_intervals=[], fine_label_unknown_intervals=[],
                   no_annotated_gene_either_strand=False)
        raw = self.seqs[contig][a:b]
        core = raw if strand == '+' else MF.rc(raw)
        tg = build_targets(win, {t['id']: t for t in local}, core)
        seg, region = tg['seg'].copy(), tg['region'].copy()
        if treatment in S_ARMS:
            n = b - a
            known = np.zeros(seg.shape[0], dtype=bool)
            for t in local:
                x, y = MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, strand)
                known[x:y] = True
            for p, q in self.background_local(contig, strand, a, b):
                x, y = MF.orient_interval(p, q, a, b, strand)
                known[x:y] = True
            known[n:] = False
            seg[~known] = IGNORE
            region[~known] = IGNORE
        return dict(seg=seg, region=region), win, core


def arm_seeds(arm, s1_rows=None, s2_rows=None):
    """Seed rows of an arm (protocol 22 §2 table): S1 / F1 -> self-training seeds, S2 -> protein seeds, S3 -> S1 ∪ S2 with S2 priority."""
    if arm in ('S1', 'F1'):
        return list(s1_rows), 0
    if arm == 'S2':
        return list(s2_rows), 0
    if arm == 'S3':
        return merge_seed_sets(s2_rows, s1_rows)
    raise ValueError(arm)


def target_labeler(species, arm, s1=None, s2=None, background=None, contigs=None):
    """Convenience constructor on the M30 label-set A data of a target species (train chromosomes by default; the validation pseudo-contig only on request)."""
    lens = W.a_lengths(species)
    if contigs is None:
        contigs = W.TARGETS[species]['train']
    lens = {c: lens[c] for c in contigs}
    rows, _ = arm_seeds(arm, s1, s2)
    return SeedLabeler(species, W.M1 / species / 'genome.fa', lens, rows, background)
