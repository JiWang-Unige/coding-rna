"""M30 labels: catalog of primary transcripts, nested labeled-gene selection, window grid, treatments U (unknown) / F (fill background).
Protocol: docs/review_20261016/13_M30_protocol.md §3.  Uses the frozen label code (`src.m28.labels.build_targets`) and the frozen ruler selection."""
import bisect, json, re
from collections import defaultdict
import numpy as np
import m30_common as M
import a1_common as K
M29 = M.register()
import common_ruler as CR                      # noqa: E402  (frozen)
import manifest as MF                          # noqa: E402  (frozen: starts, merged, orient_interval, rc, read_allowed)
import census as CE                            # noqa: E402  (frozen: overlapping_ids)
from src.m28.labels import build_targets, IGNORE   # noqa: E402

WIDTH = MF.WIDTH


class Labeler:
    def __init__(self, species, contigs, k, selection_seed=20261017):
        """contigs: contig names whose genes are catalogued; k: labeled fraction (0<k<=1) applied independently per contig list given."""
        self.species, self.contigs, self.k = species, list(contigs), k
        path = M.DATA / species
        lens = M.lengths(species)
        self.lengths = {c: lens[c] for c in self.contigs}
        self.seqs = MF.read_allowed(path / 'genome.fa', self.lengths)
        coding = CR.S.E.parse_annotation(path / 'reference.gff3', self.lengths, protein_coding_only=True)
        raw_cds = {}
        for c in self.contigs:
            with (path / 'reference.gff3').open() as f:
                raw_cds.update(CR.P.raw_features(f, c)[2])
        refs, _ = CR.select_reference(coding, raw_cds, self.seqs)
        ref_ids = {t['id'] for t in refs}
        local = CR.S.E.primary_transcripts(coding, complete_only=False)
        by_gene = {t['gene_id']: t for t in local}
        by_gene.update({t['gene_id']: t for t in refs})          # assessable primary replaces the local primary (as in frozen align_labels)
        prim = list(by_gene.values())
        overlaps = set()
        for c in self.contigs:
            overlaps |= CE.overlapping_ids([t for t in prim if t['seqid'] == c])
        cds_partial = MF.partial_cds_ids(path / 'reference.gff3', self.lengths)
        self.catalog = {}
        for t in prim:
            failures, _ = MF.chain_checks(t, self.seqs[t['seqid']])        # frozen manifest checks (ATG/stop/phase/splice motifs ...)
            supported = (t['id'] in ref_ids) and not failures and t['id'] not in cds_partial
            self.catalog[t['id']] = dict(id=t['id'], gene_id=t['gene_id'], seqid=t['seqid'], strand=t['strand'], CDS=t['CDS'],
                span=[t['CDS'][0][0], t['CDS'][-1][1]], training_primary=True, CDS_complete_supported=supported,
                overlapping_primary=t['id'] in overlaps, assessable=t['id'] in ref_ids, eligible=supported and t['id'] not in overlaps)
        # nested labeled selection among assessable genes (protocol §3): sorted by (contig, start), fixed permutation
        elig = sorted((t for t in self.catalog.values() if t['eligible']), key=lambda t: (t['seqid'], t['span'][0], t['id']))
        perm = np.random.default_rng(selection_seed).permutation(len(elig))
        n = int(np.ceil(k * len(elig)))
        self.labeled = {elig[i]['id'] for i in perm[:n]}
        self.n_eligible, self.n_labeled = len(elig), len(self.labeled)
        self.by_qs = defaultdict(list)
        for t in self.catalog.values():
            self.by_qs[t['seqid'], t['strand']].append(t)
        self.starts = {}
        for key, lst in self.by_qs.items():
            lst.sort(key=lambda t: t['span'][0])
            self.starts[key] = [t['span'][0] for t in lst]
        # max span for the left search bound
        self.max_span = {key: max((t['span'][1] - t['span'][0] for t in lst), default=0) for key, lst in self.by_qs.items()}
        genes = defaultdict(list)
        for g in CR.S.E.parse_annotation(path / 'reference.gff3', self.lengths)['genes'].values():
            genes[g['seqid']].append((g['start'], g['end']))
        self.genes = genes

    # ---- windows --------------------------------------------------------------------------------------------------
    def grid(self, contig):
        L = self.lengths[contig]
        return [(a, min(L, a + WIDTH)) for a in MF.starts(L)]

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

    def labeled_windows(self):
        """(window start, end, contig) grid cells covering each labeled gene (for gene-first sampling): returns {gene id: [(contig, strand, a, b)]}"""
        res = {}
        for tid in self.labeled:
            t = self.catalog[tid]
            c, s = t['seqid'], t['strand']
            L = self.lengths[c]
            cells = [(a, min(L, a + WIDTH)) for a in MF.starts(L)]
            starts = [a for a, _ in cells]
            i = bisect.bisect_right(starts, t['span'][0]) - 1
            out = []
            for j in range(max(0, i - 2), min(len(cells), i + 3)):
                a, b = cells[j]
                if a < t['span'][1] and b > t['span'][0]:
                    out.append((c, s, a, b))
            res[tid] = out
        return res

    # ---- targets --------------------------------------------------------------------------------------------------
    def targets(self, contig, strand, a, b, treatment):
        """treatment: 'U' (unselected-gene loci unknown) | 'F' (unselected genes filled with background) | 'FULL' (every assessable gene labeled)."""
        local = self.local(contig, strand, a, b)
        sel = {t['id'] for t in local if (t['eligible'] and (t['id'] in self.labeled or treatment == 'FULL')) or not t['eligible']}
        visible = [t for t in local if t['id'] in sel]
        win = self.window(contig, strand, a, b)
        unknown = []
        for t in visible:
            if not t['CDS_complete_supported'] or t['overlapping_primary']:
                unknown.append(MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, strand))
        win.update(overlapping_transcript_ids=[t['id'] for t in visible], positive_ids=[], chain_unknown_intervals=[],
                   fine_label_unknown_intervals=MF.merged(unknown), no_annotated_gene_either_strand=False)
        raw = self.seqs[contig][a:b]
        core = raw if strand == '+' else MF.rc(raw)
        tg = build_targets(win, {t['id']: t for t in visible}, core)
        seg, region = tg['seg'].copy(), tg['region'].copy()
        if treatment == 'U':
            for t in local:
                if t['eligible'] and t['id'] not in self.labeled:
                    x, y = MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, strand)
                    if y > x:
                        bgmask = seg[x:y] == 0                         # never erase labeled-gene roles
                        seg[x:y][bgmask] = IGNORE
                        region[x:y][region[x:y] == 0] = IGNORE
        return dict(seg=seg, region=region), win, core


def labeled_f1(gff_chains, labeler, contig, strand_free=False):
    """labeled-F1 (protocol §4): TP = exactly recovered labeled genes; ref = labeled genes; pred = predicted chains overlapping (same strand) a labeled gene."""
    lab = [t for t in labeler.catalog.values() if t['id'] in labeler.labeled and t['seqid'] == contig]
    keys = {(contig, t['strand'], tuple((a, b) for a, b, _ in t['CDS'])) for t in lab}
    by = defaultdict(list)
    for t in lab:
        by[t['strand']].append(tuple(t['span']))
    idx = {}
    for st, sp in by.items():
        sp.sort()
        starts = [x for x, _ in sp]
        pm, m = [], -1
        for _, y in sp:
            m = max(m, y); pm.append(m)
        idx[st] = (starts, pm)
    tp = pred = 0
    for k, v in gff_chains.items():
        a, b = v['exons'][0][0], v['exons'][-1][1]
        starts, pm = idx.get(v['strand'], ([], []))
        i = bisect.bisect_left(starts, b)             # labeled genes starting before the chain end
        if i > 0 and pm[i - 1] > a:                   # some labeled gene on this strand overlaps the chain
            pred += 1
            tp += int(k in keys)
    ref = len(keys)
    return dict(tp=tp, pred=pred, ref=ref, F1=100 * 2 * tp / (pred + ref) if pred + ref else 0.)
