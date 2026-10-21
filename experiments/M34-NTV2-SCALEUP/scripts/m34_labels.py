"""M34 label builder + manifest export.  Same logic as the M30 `Labeler` (itself validated against the frozen R4 targets on Arabidopsis, >99.9% agreement,
experiments/M30-HOLDOUT-ADAPT/tests) but registry-driven for all 10 M34 species, with a single GFF pass, and an export of M28-format manifests:
  data/m34_train/catalog/<species>/{transcripts.jsonl, windows_train.jsonl, windows_val.jsonl, summary.json, labeled_k25.json, labeled_k50.json}
Labels use the frozen code `src.m28.labels.build_targets` and the frozen ruler selection (`common_ruler.select_reference`).
usage: m34_labels.py build <species_key> ...        (CPU only)
"""
import bisect
import json
import sys
import time
from collections import Counter, defaultdict

import numpy as np

import m34_common as C
import a1_common as K
import common_ruler as CR                      # noqa: E402  (frozen)
import manifest as MF                          # noqa: E402  (frozen)
import census as CE                            # noqa: E402  (frozen)
from src.m28.labels import build_targets, IGNORE   # noqa: E402

WIDTH = MF.WIDTH


def group_gff_lines(path, contigs):
    by = defaultdict(list)
    keep = {'gene', 'mRNA', 'CDS'}
    with open(path) as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.split('\t', 3)
            if len(p) >= 3 and p[0] in contigs and p[2] in keep:
                by[p[0]].append(line)
    return by


class Labeler:
    def __init__(self, key, contigs, k=1.0, selection_seed=C.SELECTION_SEED):
        sp = C.species_table()[key]
        self.species, self.contigs, self.k = key, list(contigs), k
        path = sp['dir']
        self.lengths = {c: sp['lengths'][c] for c in self.contigs}
        self.seqs = MF.read_allowed(path / 'genome.fa', self.lengths)
        coding = CR.S.E.parse_annotation(path / 'reference.gff3', self.lengths, protein_coding_only=True)
        by = group_gff_lines(path / 'reference.gff3', set(self.contigs))
        raw_cds = {}
        for c in self.contigs:
            raw_cds.update(CR.P.raw_features(by[c], c)[2])
        del by
        refs, _ = CR.select_reference(coding, raw_cds, self.seqs)
        ref_ids = {t['id'] for t in refs}
        local = CR.S.E.primary_transcripts(coding, complete_only=False)
        by_gene = {t['gene_id']: t for t in local}
        by_gene.update({t['gene_id']: t for t in refs})          # assessable primary replaces the local primary (as in the frozen align_labels)
        prim = list(by_gene.values())
        overlaps = set()
        for c in self.contigs:
            overlaps |= CE.overlapping_ids([t for t in prim if t['seqid'] == c])
        cds_partial = MF.partial_cds_ids(path / 'reference.gff3', self.lengths)
        self.n_coding_genes = len(coding['genes'])
        self.catalog = {}
        for t in prim:
            failures, motifs = MF.chain_checks(t, self.seqs[t['seqid']])
            supported = (t['id'] in ref_ids) and not failures and t['id'] not in cds_partial
            self.catalog[t['id']] = dict(id=t['id'], gene_id=t['gene_id'], seqid=t['seqid'], strand=t['strand'], CDS=t['CDS'],
                span=[t['CDS'][0][0], t['CDS'][-1][1]], training_primary=True, CDS_complete_supported=supported,
                overlapping_primary=t['id'] in overlaps, assessable=t['id'] in ref_ids, CDS_assessable_primary=t['id'] in ref_ids,
                eligible=supported and t['id'] not in overlaps, eligible_positive=supported and t['id'] not in overlaps,
                sequence_structure_failures=failures, splice_motifs=motifs, CDS_partial_attribute=t['id'] in cds_partial,
                n_exons=len(t['CDS']), cds_bp=sum(b - a for a, b, _ in t['CDS']),
                split='val' if sp['val'] and t['seqid'] in sp['val'] else 'train')
        elig = [t for t in self.catalog.values() if t['eligible'] and t['seqid'] in self.contigs]
        self.labeled = C.select_labeled(elig, k, selection_seed)
        self.n_eligible, self.n_labeled = len(elig), len(self.labeled)
        self.by_qs = defaultdict(list)
        for t in self.catalog.values():
            self.by_qs[t['seqid'], t['strand']].append(t)
        self.starts, self.max_span = {}, {}
        for key2, lst in self.by_qs.items():
            lst.sort(key=lambda t: t['span'][0])
            self.starts[key2] = [t['span'][0] for t in lst]
            self.max_span[key2] = max((t['span'][1] - t['span'][0] for t in lst), default=0)
        genes = defaultdict(list)
        for g in CR.S.E.parse_annotation(path / 'reference.gff3', self.lengths)['genes'].values():
            genes[g['seqid']].append((g['start'], g['end']))
        self.genes = genes
        self.gene_index = {c: self._merge(v) for c, v in genes.items()}

    @staticmethod
    def _merge(spans):
        m = MF.merged(spans)
        return [x for x, _ in m], [y for _, y in m]

    # ---- windows --------------------------------------------------------------------------------------------------
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

    def labeled_windows(self):
        """{gene id: [(contig, strand, a, b)]} cells covering each labeled gene (same cells as M30)."""
        res = {}
        for tid in self.labeled:
            t = self.catalog[tid]
            res[tid] = [(t['seqid'], t['strand'], a, b) for a, b in C.gene_cells(t['span'], self.lengths[t['seqid']])]
        return res

    def no_gene_either_strand(self, contig, a, b):
        starts, ends = self.gene_index.get(contig, ([], []))
        i = bisect.bisect_left(ends, a + 1)             # first merged gene span with end > a
        return not (i < len(starts) and starts[i] < b)

    def window_row(self, contig, strand, a, b, split):
        """M28-manifest-format window row (label policy CDS_assessable_primary_R4)."""
        local = self.local(contig, strand, a, b)
        fully = [t for t in local if a <= t['span'][0] and t['span'][1] <= b]
        pos = [t['id'] for t in fully if t['eligible_positive'] and a + C.FLANK <= t['span'][0] and t['span'][1] <= b - C.FLANK]
        unknown, fine = [], []
        for t in local:
            iv = MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, strand)
            if not t['CDS_complete_supported'] or t['overlapping_primary']:
                fine.append(iv); unknown.append(iv)
            elif t not in fully:
                unknown.append(iv)
        raw = self.seqs[contig][a:b]
        sq = raw if strand == '+' else MF.rc(raw)
        return dict(id=f'{self.species}|{contig}|{a}|{strand}', species=self.species, split=split, seqid=contig, start=a, end=b, width=WIDTH,
                    valid_bases=b - a, strand=strand, pad_right_after_orientation=WIDTH - (b - a),
                    primary_complete_ids=[t['id'] for t in fully if t['assessable']], positive_ids=pos,
                    overlapping_transcript_ids=[t['id'] for t in local], chain_unknown_intervals=MF.merged(unknown),
                    fine_label_unknown_intervals=MF.merged(fine), no_annotated_gene_either_strand=self.no_gene_either_strand(contig, a, b),
                    motifs={'start': sq.count('ATG'), 'stop': sum(sq.count(m) for m in MF.STOP), 'donor': sq.count('GT') + sq.count('GC'), 'acceptor': sq.count('AG')},
                    label_policy='CDS_assessable_primary_R4')

    # ---- targets (identical semantics to M30) --------------------------------------------------------------------
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


def stats(lab, split_of):
    """Per-species supervision statistics over the labeler's contigs."""
    cat = [t for t in lab.catalog.values() if t['seqid'] in lab.contigs]
    out = {}
    for sp in ('train', 'val'):
        rows = [t for t in cat if split_of[t['seqid']] == sp]
        ass = [t for t in rows if t['assessable']]
        el = [t for t in rows if t['eligible']]
        bp = sum(lab.lengths[c] for c in lab.contigs if split_of[c] == sp)
        out[sp] = dict(bp=bp, primary_coding_genes=len(rows), cds_assessable_references=len(ass), eligible_positive=len(el),
            eligible_per_Mb=1e6 * len(el) / bp if bp else None,
            exon_hist={'1': sum(t['n_exons'] == 1 for t in ass), '2-5': sum(2 <= t['n_exons'] <= 5 for t in ass), '>=6': sum(t['n_exons'] >= 6 for t in ass)},
            median_exons=float(np.median([t['n_exons'] for t in ass])) if ass else None,
            spans_longer_than_window=sum(t['span'][1] - t['span'][0] > WIDTH for t in ass),
            nonstandard_splice_fraction=(sum('nonstandard_splice' in t['sequence_structure_failures'] for t in ass) / len(ass)) if ass else None,
            unsupported_reasons=dict(Counter(f for t in rows if t['assessable'] and not t['eligible'] for f in (t['sequence_structure_failures'] or ['overlap_or_partial']))))
    return out


def export(key):
    t0 = time.time()
    sp = C.species_table()[key]
    contigs = list(sp['train']) + list(sp['val'])
    split_of = {c: ('val' if c in sp['val'] else 'train') for c in contigs}
    lab = Labeler(key, contigs, 1.0)
    out = C.CATALOG / key
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'transcripts.jsonl').open('w') as f:
        for t in lab.catalog.values():
            f.write(json.dumps(t) + '\n')
    counts = {}
    for split in ('train', 'val'):
        n = 0
        with (out / f'windows_{split}.jsonl').open('w') as f:
            for c in contigs:
                if split_of[c] != split:
                    continue
                for strand in '+-':
                    for a, b in lab.grid(c):
                        f.write(json.dumps(lab.window_row(c, strand, a, b, split)) + '\n')
                        n += 1
        counts[split] = n
    train_elig = [t for t in lab.catalog.values() if t['eligible'] and split_of[t['seqid']] == 'train']
    nested = {}
    for pct in (25, 50):
        ids = sorted(C.select_labeled(train_elig, pct / 100))
        json.dump(ids, open(out / f'labeled_k{pct}.json', 'w'))
        nested[pct] = len(ids)
    summary = dict(species=key, organism=sp['organism'], accession=sp['accession'], contigs={c: lab.lengths[c] for c in contigs}, split_of=split_of,
                   n_coding_genes_parsed=lab.n_coding_genes, windows=counts, nested_labeled_counts=nested, label_policy='CDS_assessable_primary_R4 (primary-transcript catalog, as M30)',
                   selection_seed=C.SELECTION_SEED, stats=stats(lab, split_of), seconds=time.time() - t0)
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1)
    print(json.dumps({k: summary[k] for k in ('species', 'windows', 'nested_labeled_counts', 'seconds')}), json.dumps(summary['stats']['train']['eligible_positive']),
          json.dumps(summary['stats']['val']['cds_assessable_references']), flush=True)
    return summary


if __name__ == '__main__':
    K.repoint()
    if sys.argv[1] == 'build':
        for key in sys.argv[2:]:
            export(key)
