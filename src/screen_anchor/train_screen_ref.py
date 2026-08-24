"""M1-SAMEBUDGET-SCREEN-ANCHOR — train ONE random-init reference architecture under the unified
small-sample screen protocol, then predict on the held-out TEST seqids and emit per-species CDS
GFF + the test-subset genome/reference for the CDS-span evaluator.

Single seed per invocation (the sbatch loops seeds). Pure-stdlib config via argparse (no yaml dep).

Flow: load both species (FASTA+GFF) -> per-base labels -> chromosome-level split ->
train (sample_fraction of train windows) with early-stop on val macro-F1 ->
predict TEST seqids -> CDS GFF per species (outputs/<exp_id>/predictions/<species>.gff) +
test-subset genome.fa/reference.gff3 (outputs/<exp_id>/eval_subsets/<species>/).
"""
import argparse
import json
import os
import sys

import numpy as np

# allow `python -m src.screen_anchor.train_screen_ref` and direct execution
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from src.screen_anchor import data as D
    from src.screen_anchor.models import build_model
    from src.screen_anchor.gff_io import labels_to_cds_gff
    from src.screen_anchor.decoders import LinearChainCRF, LinearChainCRFVec, SemiCRF, constrained_decode
else:
    from . import data as D
    from .models import build_model
    from .gff_io import labels_to_cds_gff
    from .decoders import LinearChainCRF, LinearChainCRFVec, SemiCRF, constrained_decode


def macro_f1(conf):
    """conf: (C,C) confusion (rows=true, cols=pred). Return macro-F1 over classes."""
    f1s = []
    for c in range(conf.shape[0]):
        tp = conf[c, c]
        fp = conf[:, c].sum() - tp
        fn = conf[c, :].sum() - tp
        denom = 2 * tp + fp + fn
        f1s.append((2 * tp / denom) if denom > 0 else 0.0)
    return float(np.mean(f1s)), [float(x) for x in f1s]


def decode_batch(emissions, mask, decoder, dtype):
    """emissions (B,W,C) -> (B,W) np labels. softmax/constrained=argmax (constrained post-proc
    is applied later at the stitched-seqid level); crf=Viterbi; semicrf=segment DP."""
    if dtype in ("softmax", "constrained") or decoder is None:
        return emissions.argmax(-1).cpu().numpy()
    if dtype == "crf":
        return decoder.viterbi(emissions, mask).cpu().numpy()
    return decoder.decode(emissions, mask).cpu().numpy()  # semicrf


def evaluate(backbone, decoder, dtype, loader, device):
    """Early-stop signal = emission ARGMAX macro-F1 (fast, decoder-agnostic). The structured
    decoder (Viterbi / segment-DP) is sequential over W and too slow to run every epoch on all
    val windows; it is applied only to the FINAL test prediction (predict_seqids), where gene
    structure matters. Early-stop on emission quality is a sound proxy."""
    import torch
    backbone.eval()
    conf = np.zeros((D.NUM_CLASSES, D.NUM_CLASSES), dtype=np.int64)
    with torch.no_grad():
        for X, Y in loader:
            X = X.to(device)
            pred = backbone(X).argmax(-1).cpu().numpy()    # (B,W) emission argmax
            y = Y.numpy()
            m = y != -100
            for t, p in zip(y[m], pred[m]):
                conf[t, p] += 1
    return macro_f1(conf)


def predict_seqids(backbone, decoder, dtype, seqs, seqids, window, device,
                   min_cds_len=30, max_fill_gap=20, predict_batch=128):
    """Return {seqid: int8 predicted-class array of full length}. Decoder-aware. BATCHES windows
    across the whole test set through the backbone + decoder (the structured Viterbi loops over W
    once per BATCH, not per window) — ~predict_batch× faster than batch-1, which made the sequential
    Viterbi predict the bottleneck. For 'constrained', post-process the stitched per-seqid arrays."""
    import torch
    backbone.eval()
    if decoder is not None:
        decoder.eval()
    out = {sid: np.zeros(len(seqs[sid]), dtype=np.int8) for sid in seqids}
    # flat list of (seqid, start, end) windows across all test seqids
    idx = []
    for sid in seqids:
        L = len(seqs[sid])
        for s in range(0, L, window):
            e = min(s + window, L)
            if e - s >= 1:
                idx.append((sid, s, e))
    with torch.no_grad():
        for b0 in range(0, len(idx), predict_batch):
            chunk = idx[b0:b0 + predict_batch]
            maxw = max(e - s for _, s, e in chunk)
            X = torch.zeros(len(chunk), D.NUM_CHANNELS, maxw, dtype=torch.float32)
            M = torch.zeros(len(chunk), maxw, dtype=torch.bool)
            for i, (sid, s, e) in enumerate(chunk):
                X[i, :, :e - s] = torch.from_numpy(D.one_hot_window(seqs[sid], s, e))
                M[i, :e - s] = True
            emissions = backbone(X.to(device))             # (B, maxw, C)
            lab = decode_batch(emissions, M.to(device), decoder, dtype)  # (B, maxw)
            for i, (sid, s, e) in enumerate(chunk):
                out[sid][s:e] = lab[i, :e - s].astype(np.int8)
    if dtype == "constrained":
        out = constrained_decode(out, min_cds_len=min_cds_len, max_fill_gap=max_fill_gap)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=["tiberius_like", "helixer_like"])
    ap.add_argument("--exp-id", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--species", nargs="+", required=True,
                    help="paths to data/m1_screen/<species> dirs (each has genome.fa, reference.gff3)")
    ap.add_argument("--window", type=int, default=2048)
    ap.add_argument("--sample-fraction", type=float, default=0.3)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--patience", type=int, default=3)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-train-windows", type=int, default=0,
                    help="hard cap on train windows after sampling (0 = no cap); smoke uses small")
    ap.add_argument("--class-weighting", default="sqrt_inv", choices=["none", "inv", "sqrt_inv"],
                    help="counter class imbalance (intergenic dominates) so the model does not "
                         "collapse to one class; weights from TRAIN label frequencies.")
    ap.add_argument("--decoder", default="softmax",
                    choices=["softmax", "crf", "semicrf", "constrained"],
                    help="TA-DECODER-M3: structured decoder on the FIXED tiberius_like backbone. "
                         "softmax=per-base baseline; crf=linear-chain CRF; semicrf=semi-Markov "
                         "(max_seg_len); constrained=softmax-trained + constrained post-processing.")
    ap.add_argument("--max-seg-len", type=int, default=64, help="semicrf max segment length.")
    ap.add_argument("--crf-aux-ce", type=float, default=1.0,
                    help="aux class-weighted CE on emissions added to CRF/semiCRF NLL (counters "
                         "imbalance so emissions don't collapse, matching the weighted-CE anchor).")
    ap.add_argument("--min-cds-len", type=int, default=30, help="constrained: drop CDS runs shorter.")
    ap.add_argument("--max-fill-gap", type=int, default=20, help="constrained: fill intergenic gaps <=.")
    args = ap.parse_args()

    import torch
    from torch.utils.data import DataLoader
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    os.makedirs(args.out_dir, exist_ok=True)
    pred_dir = os.path.join(args.out_dir, "predictions")
    sub_dir = os.path.join(args.out_dir, "eval_subsets")
    os.makedirs(pred_dir, exist_ok=True)
    os.makedirs(sub_dir, exist_ok=True)

    # ---- load species ----
    species = {}  # name -> dict(seqs, labels, gff, splits)
    for sp_path in args.species:
        name = os.path.basename(sp_path.rstrip("/"))
        genome = os.path.join(sp_path, "genome.fa")
        gff = os.path.join(sp_path, "reference.gff3")
        seqs = D.read_fasta(genome)
        seq_lengths = {sid: len(s) for sid, s in seqs.items()}
        labels = D.build_labels(gff, seq_lengths)
        splits = D.assign_splits(list(seqs.keys()))
        species[name] = {"seqs": seqs, "labels": labels, "gff": gff, "splits": splits}
        cls_counts = np.bincount(
            np.concatenate([labels[s] for s in seqs]), minlength=D.NUM_CLASSES)
        print(f"[{name}] seqs={len(seqs)} splits="
              f"{ {k: sum(1 for v in splits.values() if v==k) for k in ('train','val','test')} } "
              f"label_bases(intergenic/CDS/genebody-nc)={cls_counts.tolist()}", flush=True)

    # ---- datasets (pool species per split) ----
    def pooled(split):
        seqs_all, labels_all, ids = {}, {}, []
        for name, sp in species.items():
            for sid, tag in sp["splits"].items():
                if tag == split:
                    key = f"{name}::{sid}"
                    seqs_all[key] = sp["seqs"][sid]
                    labels_all[key] = sp["labels"][sid]
                    ids.append(key)
        return seqs_all, labels_all, ids

    tr_seqs, tr_labels, tr_ids = pooled("train")
    va_seqs, va_labels, va_ids = pooled("val")
    train_ds = D.WindowDataset(tr_seqs, tr_labels, tr_ids, window=args.window,
                               sample_fraction=args.sample_fraction, seed=args.seed)
    if args.max_train_windows and len(train_ds) > args.max_train_windows:
        train_ds.index = train_ds.index[:args.max_train_windows]
    val_ds = D.WindowDataset(va_seqs, va_labels, va_ids, window=args.window,
                             sample_fraction=1.0, seed=args.seed)
    print(f"train_windows={len(train_ds)} val_windows={len(val_ds)} device={device}", flush=True)

    train_dl = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                          collate_fn=D.collate_pad)
    val_dl = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False,
                        collate_fn=D.collate_pad)

    # ---- class weights from TRAIN labels (counter intergenic dominance -> avoid collapse) ----
    cls_counts = np.zeros(D.NUM_CLASSES, dtype=np.int64)
    for sid in tr_ids:
        cls_counts += np.bincount(tr_labels[sid], minlength=D.NUM_CLASSES)
    if args.class_weighting == "none" or (cls_counts == 0).any():
        weight = None
    else:
        if args.class_weighting == "inv":
            w = 1.0 / cls_counts
        else:  # sqrt_inv
            w = 1.0 / np.sqrt(cls_counts)
        w = w / w.mean()                          # normalize to mean 1
        weight = torch.tensor(w, dtype=torch.float32, device=device)
    print(f"train_label_bases={cls_counts.tolist()} class_weights="
          f"{None if weight is None else [round(x,3) for x in weight.cpu().tolist()]}", flush=True)

    # ---- backbone + structured decoder (TA-DECODER-M3) ----
    backbone = build_model(args.model).to(device)       # TiberiusLike -> (B,W,C) emissions
    if args.decoder == "crf":
        decoder = LinearChainCRFVec(D.NUM_CLASSES).to(device)  # vectorized parallel-scan partition
    elif args.decoder == "semicrf":
        decoder = SemiCRF(D.NUM_CLASSES, max_seg_len=args.max_seg_len).to(device)
    else:
        decoder = None                                  # softmax / constrained: per-base CE
    params = list(backbone.parameters()) + (list(decoder.parameters()) if decoder else [])
    opt = torch.optim.Adam(params, lr=args.lr)
    ce = torch.nn.CrossEntropyLoss(ignore_index=-100, weight=weight)
    print(f"decoder={args.decoder} (params={sum(p.numel() for p in params)})", flush=True)

    def compute_loss(emissions, Y):
        if args.decoder in ("crf", "semicrf"):
            mask = Y != -100
            nll = decoder.nll(emissions, Y, mask)
            # aux class-weighted CE on emissions counters imbalance (match weighted-CE anchor)
            aux = ce(emissions.reshape(-1, D.NUM_CLASSES), Y.reshape(-1))
            return nll + args.crf_aux_ce * aux
        return ce(emissions.reshape(-1, D.NUM_CLASSES), Y.reshape(-1))

    best_f1, best_bb, best_dec, bad = -1.0, None, None, 0
    for ep in range(1, args.epochs + 1):
        backbone.train()
        if decoder is not None:
            decoder.train()
        tot, nb = 0.0, 0
        for X, Y in train_dl:
            X, Y = X.to(device), Y.to(device)
            opt.zero_grad()
            emissions = backbone(X)                      # (B,W,C)
            loss = compute_loss(emissions, Y)
            loss.backward()
            opt.step()
            tot += loss.item(); nb += 1
        vf1, vper = evaluate(backbone, decoder, args.decoder, val_dl, device)
        print(f"epoch {ep}: train_loss={tot/max(nb,1):.4f} val_macroF1={vf1:.4f} "
              f"per_class={ [round(x,4) for x in vper] }", flush=True)
        if vf1 > best_f1:
            best_f1, bad = vf1, 0
            best_bb = {k: v.detach().cpu().clone() for k, v in backbone.state_dict().items()}
            best_dec = (None if decoder is None
                        else {k: v.detach().cpu().clone() for k, v in decoder.state_dict().items()})
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop at epoch {ep} (best val_macroF1={best_f1:.4f})", flush=True)
                break
    if best_bb is not None:
        backbone.load_state_dict(best_bb)
        if decoder is not None and best_dec is not None:
            decoder.load_state_dict(best_dec)

    # ---- predict TEST per species -> CDS GFF + test-subset genome/ref ----
    summary = {"model": args.model, "decoder": args.decoder, "seed": args.seed,
               "best_val_macro_f1": best_f1, "device": device, "per_species_test_seqids": {}}
    for name, sp in species.items():
        test_ids = [sid for sid, tag in sp["splits"].items() if tag == "test"]
        summary["per_species_test_seqids"][name] = len(test_ids)
        if not test_ids:
            continue
        pred = predict_seqids(backbone, decoder, args.decoder, sp["seqs"], test_ids,
                              args.window, device, args.min_cds_len, args.max_fill_gap)
        n_genes = labels_to_cds_gff(pred, os.path.join(pred_dir, f"{name}.gff"),
                                    source=f"{args.model}_{args.decoder}")
        sp_sub = os.path.join(sub_dir, name)
        os.makedirs(sp_sub, exist_ok=True)
        D.write_subset_fasta(sp["seqs"], test_ids, os.path.join(sp_sub, "genome.fa"))
        D.write_subset_gff(sp["gff"], test_ids, os.path.join(sp_sub, "reference.gff3"))
        print(f"[{name}] test_seqids={len(test_ids)} predicted_genes={n_genes}", flush=True)

    with open(os.path.join(args.out_dir, "train_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("TRAIN_DONE", flush=True)


if __name__ == "__main__":
    main()
