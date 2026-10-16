"""M1-SAMEBUDGET-SCREEN-ANCHOR — data prep: genome+GFF3 -> per-base 3-class labels,
chromosome-level split, windowed one-hot dataset, and test-subset FASTA/GFF for eval.

Label scheme (per base, int8):
    0 = intergenic
    1 = CDS           (coding)
    2 = gene-body non-coding (intron / UTR)  -- "intron-ish"

Derivation (order matters; CDS wins): init all 0; for each transcript mark its
span [min(exon/CDS)..max] as 2 where currently 0; then mark CDS bases as 1.

Split is CHROMOSOME-LEVEL (by seqid) and deterministic — no window straddles the
train/test boundary and no seqid appears in two splits (leakage-safe; check_data gate).

Pure-stdlib FASTA/GFF parsing (no pyfaidx/BCBio); numpy for arrays.
"""
import gzip
import os
from collections import defaultdict

import numpy as np

# 0-based half-open internally, matching scripts/eval_gene_body_mask.py.
CLASS_INTERGENIC, CLASS_CDS, CLASS_GENEBODY_NC = 0, 1, 2
NUM_CLASSES = 3
# one-hot channels: A,C,G,T,N  (N / any non-ACGT -> channel 4, all-zero ACGT)
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
NUM_CHANNELS = 5

# --- Multi-class (M8) label scheme — SEPARATE from the 3-class screen scheme above (kept intact
# for M1-M7 reproducibility). strand-aware: CDS split by reading-frame/phase, gene-body-nc split
# into intron vs UTR, intron edges marked donor/acceptor (transition markers for the CRF). Class 0
# stays INTERGENIC so fp_penalty(genic = !=intergenic) generalizes. collapse_mc_to_3class() maps
# back to {0,1,2} so gff_io + eval are unchanged (CDS-phase->CDS; intron/UTR/donor/acceptor->genebody).
MC_INTERGENIC = 0
MC_CDS_P0, MC_CDS_P1, MC_CDS_P2 = 1, 2, 3
MC_INTRON = 4
MC_UTR = 5
MC_DONOR, MC_ACCEPTOR = 6, 7
NUM_CLASSES_MC = 8


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """seqid -> uppercase sequence string. seqid = first token after '>'."""
    seqs, cur, buf = {}, None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def _parse_attrs(attr_text):
    attrs = {}
    for part in attr_text.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v.strip().strip('"')
        elif " " in part:
            k, v = part.split(" ", 1)
            attrs[k] = v.strip().strip('"')
    return attrs


def build_labels(gff_path, seq_lengths):
    """Return {seqid: int8 array[len]} per-base labels (0/1/2). See module docstring."""
    labels = {sid: np.zeros(L, dtype=np.int8) for sid, L in seq_lengths.items()}
    # group span features by transcript (Parent of CDS/exon, or transcript_id/gene_id)
    tx_spans = defaultdict(lambda: [None, None, None])  # key -> [seqid, min0, maxEnd]
    cds_intervals = defaultdict(list)                    # seqid -> list[(s0,e)]
    with _open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            seqid, _src, feat, s, e, _sc, _st, _ph, attr = f
            if seqid not in labels:
                continue
            if feat not in ("CDS", "exon"):
                continue
            a = _parse_attrs(attr)
            gid = a.get("transcript_id") or a.get("Parent") or a.get("gene_id") or a.get("ID") or "NA"
            s0, e1 = int(s) - 1, int(e)  # 0-based half-open
            if s0 < 0 or e1 <= s0:
                continue
            key = (seqid, gid)
            span = tx_spans[key]
            span[0] = seqid
            span[1] = s0 if span[1] is None else min(span[1], s0)
            span[2] = e1 if span[2] is None else max(span[2], e1)
            if feat == "CDS":
                cds_intervals[seqid].append((s0, e1))
    # 1) paint gene-body span = 2 (only over real seqids, clipped to length)
    for (seqid, _gid), (sid2, mn, mx) in tx_spans.items():
        if sid2 is None:
            continue
        L = len(labels[seqid])
        a, b = max(0, mn), min(L, mx)
        if b > a:
            labels[seqid][a:b] = CLASS_GENEBODY_NC
    # 2) paint CDS = 1 (overrides)
    for seqid, ivs in cds_intervals.items():
        arr = labels[seqid]
        L = len(arr)
        for s0, e1 in ivs:
            a, b = max(0, s0), min(L, e1)
            if b > a:
                arr[a:b] = CLASS_CDS
    return labels


def build_labels_multiclass(gff_path, seq_lengths, splice_bp=2):
    """{seqid: int8 array} per-base MULTI-CLASS labels (M8 scheme, 8 classes). strand-aware.
    Global paint order (so CDS always wins regardless of transcript order / overlap):
      0 init intergenic -> 1 paint transcript span = INTRON -> 2 paint exon = UTR ->
      3 paint CDS = CDS_P{frame} (frame from GFF phase col, strand-aware) ->
      4 mark intron-edge splice donor/acceptor (only where still INTRON).
    Frame label = 1 + ((offset_from_5' + phase) % 3): a consistent 0/1/2 cycle (GFF phase keeps it
    continuous across a transcript's CDS segments); collapses to CDS in eval, so the exact
    biological convention only needs to be consistent + learnable (the CRF models the 3-cycle)."""
    labels = {sid: np.zeros(L, dtype=np.int8) for sid, L in seq_lengths.items()}
    # group features by transcript: spans, exon intervals, cds(start,end,phase), strand
    tx = defaultdict(lambda: {"seqid": None, "strand": "+", "exon": [], "cds": []})
    with _open(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            seqid, _src, feat, s, e, _sc, strand, phase, attr = f
            if seqid not in labels or feat not in ("CDS", "exon"):
                continue
            a = _parse_attrs(attr)
            gid = a.get("transcript_id") or a.get("Parent") or a.get("gene_id") or a.get("ID") or "NA"
            s0, e1 = int(s) - 1, int(e)
            if s0 < 0 or e1 <= s0:
                continue
            t = tx[(seqid, gid)]
            t["seqid"] = seqid
            t["strand"] = strand if strand in ("+", "-") else "+"
            if feat == "exon":
                t["exon"].append((s0, e1))
            else:
                ph = int(phase) if phase in ("0", "1", "2") else 0
                t["cds"].append((s0, e1, ph))
    spans, exons, cds_all, introns = [], [], [], []
    for (seqid, _gid), t in tx.items():
        if t["seqid"] is None:
            continue
        ivs = sorted(t["exon"]) if t["exon"] else sorted([(c[0], c[1]) for c in t["cds"]])
        if not ivs:
            continue
        mn, mx = ivs[0][0], max(e for _, e in ivs)
        spans.append((seqid, mn, mx))
        for s0, e1 in ivs:
            exons.append((seqid, s0, e1))
        for c0, c1 in zip([e for _, e in ivs[:-1]], [s for s, _ in ivs[1:]]):
            if c1 > c0:                                    # intron = gap between consecutive exons
                introns.append((seqid, c0, c1, t["strand"]))
        for s0, e1, ph in t["cds"]:
            cds_all.append((seqid, s0, e1, ph, t["strand"]))
    for seqid, mn, mx in spans:
        L = len(labels[seqid]); a, b = max(0, mn), min(L, mx)
        if b > a:
            labels[seqid][a:b] = MC_INTRON
    for seqid, s0, e1 in exons:
        L = len(labels[seqid]); a, b = max(0, s0), min(L, e1)
        if b > a:
            labels[seqid][a:b] = MC_UTR
    for seqid, s0, e1, ph, strand in cds_all:
        arr = labels[seqid]; L = len(arr)
        a, b = max(0, s0), min(L, e1)
        if b <= a:
            continue
        pos = np.arange(a, b)
        off = (pos - s0) if strand == "+" else (e1 - 1 - pos)
        frame = ((off + ph) % 3).astype(np.int8)
        arr[a:b] = (MC_CDS_P0 + frame)
    for seqid, s0, e1, strand in introns:
        arr = labels[seqid]; L = len(arr)
        d0, d1 = max(0, s0), min(L, s0 + splice_bp)           # 5' end of intron (genomic)
        a0, a1 = max(0, e1 - splice_bp), min(L, e1)            # 3' end of intron
        donor, acceptor = ((d0, d1), (a0, a1)) if strand == "+" else ((a0, a1), (d0, d1))
        if donor[1] > donor[0]:
            seg = arr[donor[0]:donor[1]]; seg[seg == MC_INTRON] = MC_DONOR
        if acceptor[1] > acceptor[0]:
            seg = arr[acceptor[0]:acceptor[1]]; seg[seg == MC_INTRON] = MC_ACCEPTOR
    return labels


def collapse_mc_to_3class(arr):
    """Map an M8 multi-class label/prediction array (0..7) -> the 3-class screen scheme
    {0 intergenic, 1 CDS, 2 gene-body-nc} so gff_io.labels_to_cds_gff + eval are unchanged."""
    out = np.zeros(len(arr), dtype=np.int8)
    a = np.asarray(arr)
    out[(a >= MC_CDS_P0) & (a <= MC_CDS_P2)] = CLASS_CDS
    out[(a == MC_INTRON) | (a == MC_UTR) | (a == MC_DONOR) | (a == MC_ACCEPTOR)] = CLASS_GENEBODY_NC
    return out


def assign_splits(seqids, val_mod=3, test_mod=4, mod=5):
    """Deterministic chromosome-level split. Sort seqids; index i: i%mod==val_mod->val,
    ==test_mod->test, else train. ~60/20/20, leakage-safe (split is by whole seqid)."""
    out = {}
    for i, sid in enumerate(sorted(seqids)):
        r = i % mod
        out[sid] = "val" if r == val_mod else ("test" if r == test_mod else "train")
    return out


def one_hot_window(seq, start, end):
    """ACGTN one-hot, shape (NUM_CHANNELS, end-start). N/other -> channel 4."""
    w = end - start
    x = np.zeros((NUM_CHANNELS, w), dtype=np.float32)
    sub = seq[start:end]
    for j, ch in enumerate(sub):
        x[BASE_TO_IDX.get(ch, 4), j] = 1.0
    return x


class WindowDataset:
    """Tiled non-overlapping windows over the given seqids. Encodes on the fly.
    Lazy import of torch so this module imports without torch (for check_data)."""

    def __init__(self, seqs, labels, seqids, window=2048, sample_fraction=1.0, seed=0):
        self.seqs, self.labels, self.window = seqs, labels, window
        self.index = []  # (seqid, start, end)
        for sid in seqids:
            L = len(seqs[sid])
            for s in range(0, L, window):
                e = min(s + window, L)
                if e - s >= 32:  # skip tiny tail windows
                    self.index.append((sid, s, e))
        if sample_fraction < 1.0 and self.index:
            rng = np.random.default_rng(seed)
            k = max(1, int(round(len(self.index) * sample_fraction)))
            sel = rng.choice(len(self.index), size=k, replace=False)
            self.index = [self.index[i] for i in sorted(sel)]

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i):
        import torch
        sid, s, e = self.index[i]
        x = one_hot_window(self.seqs[sid], s, e)            # (C, w)
        y = self.labels[sid][s:e].astype(np.int64)          # (w,)
        return torch.from_numpy(x), torch.from_numpy(y)


def collate_pad(batch):
    """Pad variable-length tail windows to the batch max; mask via label=-100 (ignored)."""
    import torch
    maxw = max(x.shape[1] for x, _ in batch)
    C = batch[0][0].shape[0]
    X = torch.zeros(len(batch), C, maxw, dtype=torch.float32)
    Y = torch.full((len(batch), maxw), -100, dtype=torch.int64)
    for i, (x, y) in enumerate(batch):
        w = x.shape[1]
        X[i, :, :w] = x
        Y[i, :w] = y
    return X, Y


def write_subset_fasta(seqs, seqids, path):
    with open(path, "w") as fh:
        for sid in seqids:
            fh.write(f">{sid}\n")
            s = seqs[sid]
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def write_subset_gff(gff_path, seqids, out_path):
    keep = set(seqids)
    with _open(gff_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.split("\t")
            if len(f) == 9 and f[0] in keep:
                out.write(line)
