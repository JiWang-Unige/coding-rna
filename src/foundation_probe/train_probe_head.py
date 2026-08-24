"""FP-SEGMENTNT-PROBE-M1 · train a LIGHT per-base head on FROZEN cached SegmentNT element
features (torch; reads the per-seqid .npz from extract_segmentnt.py). Mirrors the frozen screen
protocol of src/screen_anchor/train_screen_ref.py EXACTLY (same chromosome split via
D.assign_splits, same window tiling, sample_fraction 0.3, 8 epochs, patience 3, class-weighted
CE) so it is a SAME-BUDGET comparison vs the anchor. The ONLY change is the input: SegmentNT
present-prob over genomic elements (F channels) instead of raw one-hot DNA (the mechanism delta).

predict TEST seqids -> CDS GFF + test-subset genome/ref (same split) for eval_gene_body_mask.py.
Usage: python -m src.foundation_probe.train_probe_head --cache <dir> --species data/m1_screen/<sp> ... \
         --exp-id <id> --out-dir outputs/<id> --seed 0 [--head mlp1x1] [--sample-fraction 0.3]
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if __package__ in (None, ""):
    sys.path.insert(0, ROOT)
from src.screen_anchor import data as D                      # noqa: E402
from src.screen_anchor.gff_io import labels_to_cds_gff       # noqa: E402


_OH_LUT = np.full(256, 4, dtype=np.int64)                    # default -> channel 4 (N/other)
for _c, _i in {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}.items():
    _OH_LUT[ord(_c)] = _i


def fast_one_hot(seq, start, end):
    """Vectorized ACGTN one-hot (NUM_CHANNELS, w) — replaces D.one_hot_window's per-char Python
    loop (the DataLoader bottleneck for the fusion path). LUT + numpy indexing."""
    idx = _OH_LUT[np.frombuffer(seq[start:end].encode("ascii", "replace"), dtype=np.uint8)]
    x = np.zeros((D.NUM_CHANNELS, end - start), dtype=np.float32)
    x[idx, np.arange(end - start)] = 1.0
    return x


def fp_penalty(emissions, Y, ignore=-100):
    """Specificity-targeted intergenic-FP penalty: mean GENIC probability mass the model puts on
    bases whose TRUE label is intergenic (class 0). Minimizing it pushes predictions OFF the true
    intergenic region -> fewer intergenic false positives -> higher intergenic_specificity.
    Differentiable; complements class-weighted CE (FP-SEGNT-FPLOSS)."""
    import torch
    p = torch.softmax(emissions, dim=-1)            # (B,W,C)
    genic = 1.0 - p[..., 0]                          # P(class>0) = genic prob mass
    inter = (Y == 0)                                 # true intergenic mask (excludes pad -100)
    n = inter.sum()
    if n == 0:
        return emissions.sum() * 0.0
    return genic[inter].mean()


def macro_f1(conf):
    f1s = []
    for c in range(conf.shape[0]):
        tp = conf[c, c]; fp = conf[:, c].sum() - tp; fn = conf[c, :].sum() - tp
        den = 2 * tp + fp + fn
        f1s.append((2 * tp / den) if den > 0 else 0.0)
    return float(np.mean(f1s)), [float(x) for x in f1s]


def load_feats(cache_dir, name):
    """Return ({seqid: (L,F) float16}, features list, {seqid: split})."""
    d = np.load(os.path.join(cache_dir, f"{name}.npz"), allow_pickle=True)
    feats = {k[len("feat::"):]: d[k] for k in d.files if k.startswith("feat::")}
    splits = {sid: sp for sid, sp in zip(d["seqids"].astype(str), d["splits"].astype(str))}
    return feats, list(d["features"]), splits


class FeatureWindowDataset:
    """Tiled non-overlapping windows over given seqids; feature = cached SegmentNT (F,W) slice.
    Mirrors D.WindowDataset tiling + sampling EXACTLY (same >=32 tail rule, same rng.choice).
    If `seqs` is given (fuse-raw-dna), concat raw-DNA one-hot (5ch) ABOVE the F feature channels."""
    def __init__(self, feats, labels, seqids, window=2048, sample_fraction=1.0, seed=0, seqs=None):
        self.feats, self.labels, self.window, self.seqs = feats, labels, window, seqs
        self.index = []
        for sid in seqids:
            L = self.labels[sid].shape[0]
            for s in range(0, L, window):
                e = min(s + window, L)
                if e - s >= 32:
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
        x = self.feats[sid][s:e].astype(np.float32).T        # (F, w)
        if self.seqs is not None:
            oh = fast_one_hot(self.seqs[sid], s, e)          # (5, w) raw-DNA one-hot (vectorized)
            x = np.concatenate([x, oh], axis=0)              # (F+5, w)
        y = self.labels[sid][s:e].astype(np.int64)           # (w,)
        return torch.from_numpy(np.ascontiguousarray(x)), torch.from_numpy(y)


def collate_pad(batch):
    import torch
    maxw = max(x.shape[1] for x, _ in batch)
    C = batch[0][0].shape[0]
    X = torch.zeros(len(batch), C, maxw, dtype=torch.float32)
    Y = torch.full((len(batch), maxw), -100, dtype=torch.int64)
    for i, (x, y) in enumerate(batch):
        w = x.shape[1]
        X[i, :, :w] = x; Y[i, :w] = y
    return X, Y


class _ConvLSTMHead(__import__("torch").nn.Module):
    """SAME architecture as the anchor's TiberiusLike (Conv x2 -> biLSTM -> linear) but with
    in_channels=F (SegmentNT features) instead of 5 (one-hot DNA). This makes the probe a CLEAN
    INPUT-SIGNAL ablation: identical head/budget vs the anchor, ONLY the input differs — avoids
    confounding 'features uninformative' with 'head too weak'. Returns (B, W, C)."""
    def __init__(self, F, n_classes=None, hidden=128, conv_channels=64, lstm_layers=2, dropout=0.1):
        import torch.nn as nn
        super().__init__()
        if n_classes is None:
            n_classes = D.NUM_CLASSES
        self.conv = nn.Sequential(
            nn.Conv1d(F, conv_channels, kernel_size=9, padding=4), nn.ReLU(),
            nn.Conv1d(conv_channels, conv_channels, kernel_size=9, padding=4), nn.ReLU(),
        )
        self.lstm = nn.LSTM(conv_channels, hidden, num_layers=lstm_layers,
                            batch_first=True, bidirectional=True, dropout=dropout)
        self.head = nn.Linear(2 * hidden, n_classes)

    def forward(self, x):                 # x: (B, F, W)
        h = self.conv(x).transpose(1, 2)  # (B, W, conv_channels)
        h, _ = self.lstm(h)
        return self.head(h)               # (B, W, C)


class _Conv1dHead(__import__("torch").nn.Module):
    """Wrap a Conv1d-stack ((B,F,W)->(B,C,W)) to return (B,W,C) like the convlstm head."""
    def __init__(self, net):
        super().__init__()
        self.net = net

    def forward(self, x):
        return self.net(x).transpose(1, 2)


def build_head(name, F, n_classes=None, hidden=64):
    import torch.nn as nn
    if n_classes is None:
        n_classes = D.NUM_CLASSES
    if name == "linear":
        return _Conv1dHead(nn.Conv1d(F, n_classes, kernel_size=1))       # purest per-base ablation
    if name == "mlp1x1":
        return _Conv1dHead(nn.Sequential(                               # per-base nonlinear ablation
            nn.Conv1d(F, hidden, kernel_size=1), nn.ReLU(),
            nn.Conv1d(hidden, hidden, kernel_size=1), nn.ReLU(),
            nn.Conv1d(hidden, n_classes, kernel_size=1),
        ))
    # 'convlstm' (default): anchor-matched head -> clean input-signal ablation
    return _ConvLSTMHead(F, n_classes=n_classes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--species", nargs="+", required=True, help="data/m1_screen/<sp> dirs (genome.fa, reference.gff3)")
    ap.add_argument("--exp-id", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--head", default="convlstm", choices=["convlstm", "linear", "mlp1x1"])
    ap.add_argument("--window", type=int, default=2048)
    ap.add_argument("--sample-fraction", type=float, default=0.3)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--patience", type=int, default=3)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--class-weighting", default="sqrt_inv", choices=["none", "inv", "sqrt_inv"])
    ap.add_argument("--predict-batch", type=int, default=128)
    # --- TA-FOUNDATION-DECODER-M4: 3 orthogonal mechanisms to convert recall -> specificity ---
    ap.add_argument("--loss", default="ce", choices=["ce", "fp_aware"],
                    help="fp_aware = class-weighted CE + lambda * intergenic-FP penalty (genic prob "
                         "mass at TRUE-intergenic bases) — specificity-targeted (FP-SEGNT-FPLOSS).")
    ap.add_argument("--fp-lambda", type=float, default=1.0, help="weight of the intergenic-FP penalty.")
    ap.add_argument("--fuse-raw-dna", action="store_true",
                    help="concat raw-DNA one-hot (5ch) with frozen SegmentNT features -> head input "
                         "(FP-SEGNT-FUSION: raw-DNA specificity + foundation recall).")
    ap.add_argument("--decoder", default="none", choices=["none", "crf", "semi_crf"],
                    help="crf = LinearChainCRFVec; semi_crf = SemiCRF (segment-level) on head "
                         "emissions, trained with FP-aware aux (M8: biologically-meaningful transitions "
                         "over the strand-aware multi-class labels).")
    ap.add_argument("--label-scheme", default="3class", choices=["3class", "multiclass"],
                    help="3class = {intergenic,CDS,genebody-nc} (M1-M7). multiclass = M8 8-class "
                         "strand-aware {intergenic,CDS-phase0/1/2,intron,UTR,donor,acceptor}; collapses "
                         "to 3class for gff/eval (so the CDS/intergenic metrics are unchanged).")
    ap.add_argument("--semi-crf-max-seg", type=int, default=64, help="SemiCRF bounded segment length.")
    ap.add_argument("--crf-aux-ce", type=float, default=1.0, help="aux CE weight added to CRF NLL.")
    # --- TA-COHERENCE-FIX-M5: deterministic post-processing to de-fragment predictions ---
    ap.add_argument("--postproc", default="none", choices=["none", "constrained"],
                    help="constrained = src/screen_anchor/decoders.constrained_decode on the predicted "
                         "per-seqid arrays before GFF (merge small intergenic gaps / drop tiny CDS) to "
                         "cut gene_count_ratio (fragmentation) WITHOUT a learned CRF's instability.")
    ap.add_argument("--min-cds-len", type=int, default=30, help="constrained: drop CDS runs shorter.")
    ap.add_argument("--max-fill-gap", type=int, default=20, help="constrained: fill intergenic gaps <=.")
    ap.add_argument("--save-raw-pred", action="store_true",
                    help="M6: save the RAW (pre-constrained) per-seqid prediction arrays for VAL + TEST "
                         "to outputs/<exp>/raw_pred/{split}_{species}.npz, for an offline constrained-param "
                         "sweep chosen on VAL (no test leakage). Also writes val eval_subsets.")
    args = ap.parse_args()

    import torch
    from torch.utils.data import DataLoader
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    # label scheme: 3class (M1-M7) vs multiclass (M8). NCLS sizes head/weights/decoder/conf; the
    # multi-class predictions COLLAPSE to 3class for gff_io + eval (CDS-phase->CDS, intron/UTR/splice
    # ->genebody) so the CDS/intergenic metrics are identical-ruler vs prior rounds.
    MULTI = args.label_scheme == "multiclass"
    NCLS = D.NUM_CLASSES_MC if MULTI else D.NUM_CLASSES
    build_lbl = D.build_labels_multiclass if MULTI else D.build_labels
    collapse = (lambda a: D.collapse_mc_to_3class(a)) if MULTI else (lambda a: a)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(args.out_dir, exist_ok=True)
    pred_dir = os.path.join(args.out_dir, "predictions"); os.makedirs(pred_dir, exist_ok=True)
    sub_dir = os.path.join(args.out_dir, "eval_subsets"); os.makedirs(sub_dir, exist_ok=True)

    # ---- load species: genome (for eval_subsets), labels, split, cached features ----
    species, F, features = {}, None, None
    for sp_path in args.species:
        name = os.path.basename(sp_path.rstrip("/"))
        seqs = D.read_fasta(os.path.join(sp_path, "genome.fa"))
        seq_lengths = {sid: len(s) for sid, s in seqs.items()}
        labels = build_lbl(os.path.join(sp_path, "reference.gff3"), seq_lengths)
        feats, feat_names, splits = load_feats(args.cache, name)
        if features is None:
            features, F = feat_names, len(feat_names)
        assert feat_names == features, "feature order mismatch across species caches"
        # sanity: cached feat length must equal seqid length (alignment guard)
        for sid in feats:
            assert feats[sid].shape[0] == seq_lengths[sid], (
                f"{name}::{sid} feat L={feats[sid].shape[0]} != genome L={seq_lengths[sid]}")
        species[name] = {"sp_path": sp_path, "seqs": seqs, "labels": labels,
                         "feats": feats, "splits": splits}
        print(f"[{name}] seqs={len(seqs)} F={F} splits="
              f"{ {k: sum(1 for v in splits.values() if v==k) for k in ('train','val','test')} }", flush=True)

    # ---- pool species per split (key 'name::seqid', mirrors train_screen_ref) ----
    def pooled(split):
        feats_all, labels_all, seqs_all, ids = {}, {}, {}, []
        for name, sp in species.items():
            for sid, tag in sp["splits"].items():
                if tag == split and sid in sp["feats"]:
                    key = f"{name}::{sid}"
                    feats_all[key] = sp["feats"][sid]
                    labels_all[key] = sp["labels"][sid]
                    seqs_all[key] = sp["seqs"][sid]
                    ids.append(key)
        return feats_all, labels_all, seqs_all, ids

    tr_f, tr_l, tr_s, tr_ids = pooled("train")
    va_f, va_l, va_s, va_ids = pooled("val")
    fuse_seqs_tr = tr_s if args.fuse_raw_dna else None
    fuse_seqs_va = va_s if args.fuse_raw_dna else None
    train_ds = FeatureWindowDataset(tr_f, tr_l, tr_ids, window=args.window,
                                    sample_fraction=args.sample_fraction, seed=args.seed, seqs=fuse_seqs_tr)
    val_ds = FeatureWindowDataset(va_f, va_l, va_ids, window=args.window, sample_fraction=1.0,
                                  seed=args.seed, seqs=fuse_seqs_va)
    F_in = F + (D.NUM_CHANNELS if args.fuse_raw_dna else 0)
    print(f"train_windows={len(train_ds)} val_windows={len(val_ds)} device={device} "
          f"F_in={F_in} (fuse_raw_dna={args.fuse_raw_dna})", flush=True)
    train_dl = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, collate_fn=collate_pad)
    val_dl = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, collate_fn=collate_pad)

    # ---- class weights from TRAIN labels (counter intergenic dominance) ----
    cls = np.zeros(NCLS, dtype=np.int64)
    for sid in tr_ids:
        cls += np.bincount(tr_l[sid], minlength=NCLS)
    if args.class_weighting == "none" or (cls == 0).any():
        weight = None
    else:
        ww = (1.0 / cls) if args.class_weighting == "inv" else (1.0 / np.sqrt(cls))
        ww = ww / ww.mean()
        weight = torch.tensor(ww, dtype=torch.float32, device=device)
    print(f"train_label_bases={cls.tolist()} class_weights="
          f"{None if weight is None else [round(x,3) for x in weight.cpu().tolist()]}", flush=True)

    net = build_head(args.head, F_in, n_classes=NCLS).to(device)
    decoder = None
    if args.decoder == "crf":
        from src.screen_anchor.decoders import LinearChainCRFVec
        decoder = LinearChainCRFVec(NCLS).to(device)
    elif args.decoder == "semi_crf":
        from src.screen_anchor.decoders import SemiCRF
        decoder = SemiCRF(NCLS, max_seg_len=args.semi_crf_max_seg).to(device)
    params = list(net.parameters()) + (list(decoder.parameters()) if decoder else [])
    opt = torch.optim.Adam(params, lr=args.lr)
    ce = torch.nn.CrossEntropyLoss(ignore_index=-100, weight=weight)
    print(f"head={args.head} decoder={args.decoder} loss={args.loss}(lambda={args.fp_lambda}) "
          f"params={sum(p.numel() for p in params)}", flush=True)

    def fwd(X):                                              # (B,F,W) -> (B,W,C) emissions (head returns BWC)
        return net(X)

    def compute_loss(emissions, Y):
        if decoder is not None:                              # CRF: NLL + aux class-weighted CE
            mask = Y != -100
            base = decoder.nll(emissions, Y, mask) + args.crf_aux_ce * ce(
                emissions.reshape(-1, NCLS), Y.reshape(-1))
        else:
            base = ce(emissions.reshape(-1, NCLS), Y.reshape(-1))
        if args.loss == "fp_aware":
            base = base + args.fp_lambda * fp_penalty(emissions, Y)
        return base

    def evaluate():                                          # early-stop on emission ARGMAX macro-F1
        net.eval()                                           # (decoder-agnostic, fast; same as anchor)
        # multiclass: COLLAPSE pred+target to 3class so early-stop tracks the REAL objective
        # (CDS/intergenic/genebody gbF1), not the rare phase/splice classes that drag 8-class macro.
        ncf = D.NUM_CLASSES if MULTI else NCLS
        conf = np.zeros((ncf, ncf), dtype=np.int64)
        with torch.no_grad():
            for X, Y in val_dl:
                pred = fwd(X.to(device)).argmax(-1).cpu().numpy()
                y = Y.numpy(); m = y != -100
                yt, pt = y[m], pred[m]
                if MULTI:
                    yt, pt = collapse(yt), collapse(pt)
                for t, p in zip(yt, pt):
                    conf[t, p] += 1
        return macro_f1(conf)

    best_f1, best_state, best_dec, bad = -1.0, None, None, 0
    for ep in range(1, args.epochs + 1):
        net.train()
        if decoder is not None:
            decoder.train()
        tot, nb = 0.0, 0
        for X, Y in train_dl:
            X, Y = X.to(device), Y.to(device)
            opt.zero_grad()
            loss = compute_loss(fwd(X), Y)
            loss.backward(); opt.step()
            tot += loss.item(); nb += 1
        vf1, vper = evaluate()
        print(f"epoch {ep}: train_loss={tot/max(nb,1):.4f} val_macroF1={vf1:.4f} "
              f"per_class={[round(x,4) for x in vper]}", flush=True)
        if vf1 > best_f1:
            best_f1, bad = vf1, 0
            best_state = {kk: vv.detach().cpu().clone() for kk, vv in net.state_dict().items()}
            best_dec = (None if decoder is None
                        else {kk: vv.detach().cpu().clone() for kk, vv in decoder.state_dict().items()})
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop at epoch {ep} (best val_macroF1={best_f1:.4f})", flush=True)
                break
    if best_state is not None:
        net.load_state_dict(best_state)
        if decoder is not None and best_dec is not None:
            decoder.load_state_dict(best_dec)

    # ---- predict TEST per species -> CDS GFF + test-subset genome/ref ----
    net.eval()
    summary = {"head": args.head, "decoder": args.decoder, "loss": args.loss,
               "fp_lambda": args.fp_lambda, "fuse_raw_dna": args.fuse_raw_dna,
               "seed": args.seed, "best_val_macro_f1": best_f1, "device": device,
               "per_species_test_seqids": {}}
    _src = f"fp_probe_{args.head}_{args.decoder}_{args.loss}_pp-{args.postproc}"
    if decoder is not None:
        decoder.eval()
    raw_dir = os.path.join(args.out_dir, "raw_pred")
    if args.save_raw_pred:
        os.makedirs(raw_dir, exist_ok=True)

    def predict_raw(sp, seqids):
        """RAW (pre-constrained) per-seqid class arrays via batched predict over the given seqids."""
        out = {sid: np.zeros(len(sp["seqs"][sid]), dtype=np.int8) for sid in seqids}
        idx = [(sid, s, min(s + args.window, len(sp["seqs"][sid])))
               for sid in seqids for s in range(0, len(sp["seqs"][sid]), args.window)
               if min(s + args.window, len(sp["seqs"][sid])) - s >= 1]
        with torch.no_grad():
            for b0 in range(0, len(idx), args.predict_batch):
                chunk = idx[b0:b0 + args.predict_batch]
                maxw = max(e - s for _, s, e in chunk)
                X = torch.zeros(len(chunk), F_in, maxw, dtype=torch.float32)
                M = torch.zeros(len(chunk), maxw, dtype=torch.bool)
                for i, (sid, s, e) in enumerate(chunk):
                    feat = sp["feats"][sid][s:e].astype(np.float32).T
                    if args.fuse_raw_dna:
                        feat = np.concatenate([feat, fast_one_hot(sp["seqs"][sid], s, e)], axis=0)
                    X[i, :, :e - s] = torch.from_numpy(np.ascontiguousarray(feat)); M[i, :e - s] = True
                emissions = fwd(X.to(device))
                if decoder is None:
                    lab = emissions.argmax(-1).cpu().numpy()
                elif args.decoder == "semi_crf":
                    lab = decoder.decode(emissions, M.to(device)).cpu().numpy()
                else:
                    lab = decoder.viterbi(emissions, M.to(device)).cpu().numpy()
                for i, (sid, s, e) in enumerate(chunk):
                    out[sid][s:e] = lab[i, :e - s].astype(np.int8)
        # multiclass -> collapse to 3class {intergenic,CDS,genebody} so constrained_decode +
        # labels_to_cds_gff + eval are identical-ruler vs prior rounds (raw_pred saved as 3class too).
        if MULTI:
            out = {sid: collapse(arr) for sid, arr in out.items()}
        return out

    def write_subsets(sp, ids, root):
        d = os.path.join(root, name); os.makedirs(d, exist_ok=True)
        D.write_subset_fasta(sp["seqs"], ids, os.path.join(d, "genome.fa"))
        D.write_subset_gff(os.path.join(sp["sp_path"], "reference.gff3"), ids, os.path.join(d, "reference.gff3"))

    for name, sp in species.items():
        test_ids = [sid for sid, tag in sp["splits"].items() if tag == "test" and sid in sp["feats"]]
        summary["per_species_test_seqids"][name] = len(test_ids)
        if not test_ids:
            continue
        pred = predict_raw(sp, test_ids)                              # RAW test predictions
        if args.save_raw_pred:                                       # M6: cache raw for offline VAL-chosen sweep
            np.savez_compressed(os.path.join(raw_dir, f"test_{name}.npz"), **{sid: pred[sid] for sid in test_ids})
            val_ids = [sid for sid, tag in sp["splits"].items() if tag == "val" and sid in sp["feats"]]
            if val_ids:
                vpred = predict_raw(sp, val_ids)
                np.savez_compressed(os.path.join(raw_dir, f"val_{name}.npz"), **{sid: vpred[sid] for sid in val_ids})
                write_subsets(sp, val_ids, os.path.join(args.out_dir, "val_eval_subsets"))
        if args.postproc == "constrained":                            # de-fragment (M5): merge gaps / drop tiny CDS
            from src.screen_anchor.decoders import constrained_decode
            pred = constrained_decode(pred, min_cds_len=args.min_cds_len, max_fill_gap=args.max_fill_gap)
        n_genes = labels_to_cds_gff(pred, os.path.join(pred_dir, f"{name}.gff"), source=_src)
        write_subsets(sp, test_ids, sub_dir)
        print(f"[{name}] test_seqids={len(test_ids)} predicted_genes={n_genes}", flush=True)

    with open(os.path.join(args.out_dir, "train_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("TRAIN_DONE", flush=True)


if __name__ == "__main__":
    main()
