"""TB-UNFREEZE-BACKBONE-M9 · train a per-base head on a (partially-UNFROZEN) NT-v2-500m backbone.
Runs in the `generanno` env (transformers 4.49; coding-rna's 5.11 breaks the custom modeling_esm —
find_pruneable_heads_and_indices removed in 5.x; see docs/10 2026-06-12). Tests whether unfreezing
the FOUNDATION backbone lifts gene_body_F1 past the frozen-feature ceiling (M8: frozen capped ~0.74
vs ANNEVO end-to-end 0.8976). Same chromosome split / window / sample protocol as the frozen probe
(src/foundation_probe/train_probe_head.py) so it is a same-budget comparison.

  backbone = AutoModelForMaskedLM(NT-v2-500m, trust_remote_code).esm  (EsmModel, 29 layers, 1024)
  freeze all -> unfreeze top --unfreeze-layers (0 = frozen-backbone CONTROL = backbone-only self-train)
  per 2046bp window (6-divisible -> exact 6-mer alignment): tokenize -> backbone -> strip special
  tokens -> repeat_interleave x6 to per-base (1024,W) -> reuse the anchor-matched 3c FP-aware head
  -> per-base 3-class + intergenic-FP penalty -> constrained post-proc -> CDS GFF -> eval (same ruler).

Usage: python -m src.foundation_probe.train_unfreeze_backbone --species data/m1_screen/<sp> ... \
  --exp-id <id> --out-dir outputs/<id> --unfreeze-layers {0,2,4} --seed 0 [--loss fp_aware ...]
"""
import argparse, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if __package__ in (None, ""):
    sys.path.insert(0, ROOT)
from src.screen_anchor import data as D                            # noqa: E402
from src.screen_anchor.gff_io import labels_to_cds_gff             # noqa: E402
from src.foundation_probe.train_probe_head import (                # noqa: E402  reuse head + loss
    _ConvLSTMHead, fp_penalty, macro_f1)

MODEL = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
KMER = 6                                                           # NT 6-mer non-overlapping tokenizer
import re as _re                                                   # noqa: E402
_NONACGT = _re.compile(r"[^ACGT]")


def _clean(seq):
    """Uppercase + map any non-ACGT (N, softmask, IUPAC) -> A so the 6-mer tokenizer yields a
    UNIFORM token count per W-bp window (N-containing 6-mers otherwise tokenize variably -> breaks
    the fixed x6 token<->base alignment). Mirrors the N->A handling in extract_segmentnt.py."""
    return _NONACGT.sub("A", seq.upper())


class WindowTokDataset:
    """Full-W windows only (uniform token count -> no token padding in a batch). Yields (input_ids,
    per-base labels). W must be a multiple of KMER for exact token<->base alignment."""
    def __init__(self, seqs, labels, seqids, tokenizer, window, sample_fraction=1.0, seed=0, max_windows=None):
        import torch  # noqa
        self.seqs, self.labels, self.tok, self.window = seqs, labels, tokenizer, window
        self.index = []
        for sid in seqids:
            L = len(seqs[sid])
            for s in range(0, L, window):
                if s + window <= L:                                # FULL windows only
                    self.index.append((sid, s, s + window))
        if sample_fraction < 1.0 and self.index:
            rng = np.random.default_rng(seed)
            k = max(1, int(round(len(self.index) * sample_fraction)))
            sel = rng.choice(len(self.index), size=k, replace=False)
            self.index = [self.index[i] for i in sorted(sel)]
        if max_windows is not None and max_windows > 0 and len(self.index) > max_windows:
            self.index = self.index[:max_windows]

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i):
        import torch
        sid, s, e = self.index[i]
        ids = self.tok(_clean(self.seqs[sid][s:e]), return_tensors="pt")["input_ids"][0]
        y = self.labels[sid][s:e].astype(np.int64)
        return ids, torch.from_numpy(y)


def collate(batch):
    import torch
    ids = torch.stack([b[0] for b in batch])      # (B,T) uniform
    Y = torch.stack([b[1] for b in batch])         # (B,W)
    return ids, Y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--train-species-names", nargs="*", default=None,
                    help="Optional species-name allowlist for TRAIN seqids. Names are data directory basenames. "
                         "Default: all provided species. Used by M12A fixed-model cross-species protocol.")
    ap.add_argument("--val-species-names", nargs="*", default=None,
                    help="Optional species-name allowlist for VAL seqids/raw-score calibration. Default: all provided species.")
    ap.add_argument("--test-species-names", nargs="*", default=None,
                    help="Optional species-name allowlist for TEST prediction/evaluation subsets. Default: all provided species.")
    ap.add_argument("--exp-id", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--unfreeze-layers", type=int, default=2, help="0=frozen-backbone control; >0 unfreeze top N esm layers")
    ap.add_argument("--window", type=int, default=2046, help="must be multiple of 6 (kmer) for exact alignment")
    ap.add_argument("--sample-fraction", type=float, default=0.3)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--patience", type=int, default=3)
    ap.add_argument("--batch-size", type=int, default=4, help="small: backbone is heavy")
    ap.add_argument("--lr", type=float, default=1e-3, help="head LR")
    ap.add_argument("--backbone-lr", type=float, default=1e-5, help="low LR for unfrozen backbone layers")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--class-weighting", default="sqrt_inv", choices=["none", "inv", "sqrt_inv"])
    ap.add_argument("--loss", default="fp_aware", choices=["ce", "fp_aware"])
    ap.add_argument("--fp-lambda", type=float, default=1.0)
    ap.add_argument("--postproc", default="constrained", choices=["none", "constrained"])
    ap.add_argument("--min-cds-len", type=int, default=60)
    ap.add_argument("--max-fill-gap", type=int, default=20)
    ap.add_argument("--predict-batch", type=int, default=8)
    ap.add_argument("--limit-train-windows", type=int, default=None,
                    help="Optional cap after sample_fraction for large multi-species diagnostics.")
    ap.add_argument("--limit-val-windows", type=int, default=None,
                    help="Optional validation-window cap for large multi-species diagnostics.")
    ap.add_argument("--save-raw-scores", action="store_true",
                    help="Save pre-postproc per-base emissions for VAL and TEST to raw_scores/*.npz. "
                         "Used by M11 validation-only decode/FPR calibration; does not choose params on test.")
    args = ap.parse_args()
    assert args.window % KMER == 0, f"--window {args.window} must be a multiple of {KMER}"

    import torch
    from torch.utils.data import DataLoader
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(args.out_dir, exist_ok=True)
    pred_dir = os.path.join(args.out_dir, "predictions"); os.makedirs(pred_dir, exist_ok=True)
    sub_dir = os.path.join(args.out_dir, "eval_subsets"); os.makedirs(sub_dir, exist_ok=True)
    NCLS = D.NUM_CLASSES

    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
    backbone = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True).esm.to(device)
    try:
        backbone.gradient_checkpointing_enable()                   # optional VRAM saver
        print("gradient_checkpointing enabled", flush=True)
    except (ValueError, AttributeError) as _e:
        print(f"gradient_checkpointing NOT supported ({_e}); continuing without (top-N unfreeze "
              f"+ bf16 + small batch fits 24GB)", flush=True)
    hidden_dim = backbone.config.hidden_size
    # freeze all, then unfreeze top N encoder layers
    for p in backbone.parameters():
        p.requires_grad_(False)
    enc = backbone.encoder.layer
    n_unfreeze = max(0, min(args.unfreeze_layers, len(enc)))
    for layer in enc[len(enc) - n_unfreeze:]:
        for p in layer.parameters():
            p.requires_grad_(True)
    n_train_bb = sum(p.numel() for p in backbone.parameters() if p.requires_grad)
    print(f"backbone {type(backbone).__name__} layers={len(enc)} unfreeze_top={n_unfreeze} "
          f"trainable_bb_params={n_train_bb} hidden={hidden_dim}", flush=True)

    # special-token layout (probe once; all full-W windows share it -> uniform content slice)
    probe = tok("A" * args.window, return_special_tokens_mask=True)
    content_pos = [i for i, m in enumerate(probe["special_tokens_mask"]) if m == 0]
    assert len(content_pos) * KMER == args.window, (
        f"content tokens {len(content_pos)} x{KMER} != window {args.window}; check tokenizer")
    content_pos_t = torch.tensor(content_pos, device=device)
    print(f"token alignment: {len(probe['input_ids'])} tokens, {len(content_pos)} content x{KMER} = {args.window} bp", flush=True)

    # ---- load species (genome/labels/split), pool per split ----
    species = {}
    for sp_path in args.species:
        name = os.path.basename(sp_path.rstrip("/"))
        seqs = D.read_fasta(os.path.join(sp_path, "genome.fa"))
        seq_lengths = {sid: len(s) for sid, s in seqs.items()}
        labels = D.build_labels(os.path.join(sp_path, "reference.gff3"), seq_lengths)
        splits = D.assign_splits(list(seqs.keys()))
        species[name] = {"sp_path": sp_path, "seqs": seqs, "labels": labels, "splits": splits}
        print(f"[{name}] seqs={len(seqs)} splits="
              f"{ {k: sum(1 for v in splits.values() if v==k) for k in ('train','val','test')} }", flush=True)

    def _allowlist(values, known, label):
        if values is None:
            return set(known)
        allowed = set(values)
        unknown = sorted(allowed - set(known))
        if unknown:
            raise ValueError(f"unknown {label} species names: {unknown}; known={sorted(known)}")
        if not allowed:
            raise ValueError(f"{label} species allowlist is empty")
        return allowed

    known_species = set(species)
    train_species = _allowlist(args.train_species_names, known_species, "train")
    val_species = _allowlist(args.val_species_names, known_species, "val")
    test_species = _allowlist(args.test_species_names, known_species, "test")
    print(f"split_species_policy train={sorted(train_species)} val={sorted(val_species)} "
          f"test={sorted(test_species)}", flush=True)

    def pooled(split, allowed_species):
        s_, l_, ids = {}, {}, []
        for name, sp in species.items():
            if name not in allowed_species:
                continue
            for sid, tag in sp["splits"].items():
                if tag == split:
                    key = f"{name}::{sid}"
                    s_[key] = sp["seqs"][sid]; l_[key] = sp["labels"][sid]; ids.append(key)
        return s_, l_, ids

    tr_s, tr_l, tr_ids = pooled("train", train_species)
    va_s, va_l, va_ids = pooled("val", val_species)
    if not tr_ids:
        raise ValueError(f"no TRAIN seqids selected for train_species={sorted(train_species)}")
    if not va_ids:
        raise ValueError(f"no VAL seqids selected for val_species={sorted(val_species)}")
    train_ds = WindowTokDataset(
        tr_s, tr_l, tr_ids, tok, args.window, args.sample_fraction, args.seed,
        max_windows=args.limit_train_windows)
    val_ds = WindowTokDataset(
        va_s, va_l, va_ids, tok, args.window, 1.0, args.seed,
        max_windows=args.limit_val_windows)
    print(f"train_windows={len(train_ds)} val_windows={len(val_ds)} device={device}", flush=True)
    train_dl = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, collate_fn=collate)
    val_dl = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, collate_fn=collate)

    # class weights from TRAIN
    cls = np.zeros(NCLS, dtype=np.int64)
    for sid in tr_ids:
        cls += np.bincount(tr_l[sid], minlength=NCLS)
    if args.class_weighting == "none" or (cls == 0).any():
        weight = None
    else:
        ww = (1.0 / cls) if args.class_weighting == "inv" else (1.0 / np.sqrt(cls))
        weight = torch.tensor((ww / ww.mean()), dtype=torch.float32, device=device)

    head = _ConvLSTMHead(hidden_dim, n_classes=NCLS).to(device)
    params = [{"params": head.parameters(), "lr": args.lr}]
    if n_unfreeze > 0:
        params.append({"params": [p for p in backbone.parameters() if p.requires_grad], "lr": args.backbone_lr})
    opt = torch.optim.Adam(params)
    ce = torch.nn.CrossEntropyLoss(ignore_index=-100, weight=weight)
    print(f"head=convlstm unfreeze={n_unfreeze} loss={args.loss} params(head)={sum(p.numel() for p in head.parameters())}", flush=True)

    def fwd(ids):                                                  # ids (B,T) -> (B,W,NCLS)
        out = backbone(input_ids=ids)
        h = out.last_hidden_state                                  # (B,T,hidden)
        content = h.index_select(1, content_pos_t)                 # (B, W/6, hidden)
        up = content.repeat_interleave(KMER, dim=1)                # (B, W, hidden)
        return head(up.transpose(1, 2))                            # head wants (B,F,W) -> (B,W,NCLS)

    def compute_loss(em, Y):
        base = ce(em.reshape(-1, NCLS), Y.reshape(-1))
        if args.loss == "fp_aware":
            base = base + args.fp_lambda * fp_penalty(em, Y)
        return base

    def evaluate():
        head.eval(); backbone.eval()
        conf = np.zeros((NCLS, NCLS), dtype=np.int64)
        with torch.no_grad():
            for ids, Y in val_dl:
                pred = fwd(ids.to(device)).argmax(-1).cpu().numpy()
                y = Y.numpy(); m = y != -100
                for t, p in zip(y[m], pred[m]):
                    conf[t, p] += 1
        return macro_f1(conf)

    use_amp = device == "cuda"
    scaler = torch.cuda.amp.GradScaler(enabled=False)              # bf16 needs no scaler
    best_f1, best_h, best_b, bad = -1.0, None, None, 0
    for ep in range(1, args.epochs + 1):
        head.train(); backbone.train()
        tot, nb = 0.0, 0
        for ids, Y in train_dl:
            ids, Y = ids.to(device), Y.to(device)
            opt.zero_grad()
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                loss = compute_loss(fwd(ids), Y)
            loss.backward(); opt.step()
            tot += loss.item(); nb += 1
        vf1, vper = evaluate()
        print(f"epoch {ep}: train_loss={tot/max(nb,1):.4f} val_macroF1={vf1:.4f} per_class={[round(x,4) for x in vper]}", flush=True)
        if vf1 > best_f1:
            best_f1, bad = vf1, 0
            best_h = {k: v.detach().cpu().clone() for k, v in head.state_dict().items()}
            if n_unfreeze > 0:
                best_b = {k: v.detach().cpu().clone() for k, v in backbone.state_dict().items()}
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop ep {ep} (best {best_f1:.4f})", flush=True); break
    if best_h is not None:
        head.load_state_dict(best_h)
        if best_b is not None:
            backbone.load_state_dict(best_b)

    # ---- predict TEST -> 3class arrays -> constrained -> CDS GFF + subsets ----
    head.eval(); backbone.eval()
    _src = f"unfreeze{n_unfreeze}_{args.loss}_pp-{args.postproc}"
    raw_score_dir = os.path.join(args.out_dir, "raw_scores")
    if args.save_raw_scores:
        os.makedirs(raw_score_dir, exist_ok=True)

    def predict_scores(sp, seqids):
        out = {sid: np.zeros((len(sp["seqs"][sid]), NCLS), dtype=np.float16) for sid in seqids}
        idx = [(sid, s, s + args.window) for sid in seqids
               for s in range(0, len(sp["seqs"][sid]), args.window) if s + args.window <= len(sp["seqs"][sid])]
        print(f"predict_scores seqids={len(seqids)} windows={len(idx)} batch={args.predict_batch}", flush=True)
        with torch.no_grad():
            for b0 in range(0, len(idx), args.predict_batch):
                chunk = idx[b0:b0 + args.predict_batch]
                ids = torch.stack([tok(_clean(sp["seqs"][sid][s:e]), return_tensors="pt")["input_ids"][0]
                                   for sid, s, e in chunk]).to(device)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                    score = fwd(ids).float().cpu().numpy()
                for i, (sid, s, e) in enumerate(chunk):
                    out[sid][s:e, :] = score[i].astype(np.float16)
                if b0 == 0 or (b0 // args.predict_batch + 1) % 500 == 0 or b0 + args.predict_batch >= len(idx):
                    done = min(b0 + args.predict_batch, len(idx))
                    print(f"predict_scores progress {done}/{len(idx)}", flush=True)
        return out

    def labels_from_scores(score_by_seqid, intergenic_bias=0.0):
        pred = {}
        for sid, score in score_by_seqid.items():
            adjusted = score.astype(np.float32, copy=True)
            adjusted[:, 0] += intergenic_bias
            pred[sid] = adjusted.argmax(axis=-1).astype(np.int8)
        return pred

    def save_scores(split, name, score_by_seqid):
        if not args.save_raw_scores:
            return
        payload = {"seqids": np.array(sorted(score_by_seqid), dtype=str)}
        payload.update({f"score::{sid}": score_by_seqid[sid] for sid in sorted(score_by_seqid)})
        np.savez(os.path.join(raw_score_dir, f"{split}_{name}.npz"), **payload)

    def write_subsets(sp, ids, root, name):
        d = os.path.join(root, name); os.makedirs(d, exist_ok=True)
        D.write_subset_fasta(sp["seqs"], ids, os.path.join(d, "genome.fa"))
        D.write_subset_gff(os.path.join(sp["sp_path"], "reference.gff3"), ids, os.path.join(d, "reference.gff3"))

    summary = {"unfreeze_layers": n_unfreeze, "loss": args.loss, "seed": args.seed,
               "best_val_macro_f1": best_f1, "device": device, "window": args.window,
               "raw_scores_saved": bool(args.save_raw_scores),
               "train_species_names": sorted(train_species),
               "val_species_names": sorted(val_species),
               "test_species_names": sorted(test_species)}
    for name, sp in species.items():
        if args.save_raw_scores and name in val_species:
            val_ids = [sid for sid, tag in sp["splits"].items() if tag == "val"]
            if val_ids:
                val_scores = predict_scores(sp, val_ids)
                save_scores("val", name, val_scores)
                write_subsets(sp, val_ids, sub_dir, f"val_{name}")
        if name not in test_species:
            continue
        test_ids = [sid for sid, tag in sp["splits"].items() if tag == "test"]
        if not test_ids:
            continue
        test_scores = predict_scores(sp, test_ids)
        save_scores("test", name, test_scores)
        pred = labels_from_scores(test_scores)
        if args.postproc == "constrained":
            from src.screen_anchor.decoders import constrained_decode
            pred = constrained_decode(pred, min_cds_len=args.min_cds_len, max_fill_gap=args.max_fill_gap)
        ng = labels_to_cds_gff(pred, os.path.join(pred_dir, f"{name}.gff"), source=_src)
        write_subsets(sp, test_ids, sub_dir, name)
        print(f"[{name}] test_seqids={len(test_ids)} predicted_genes={ng}", flush=True)
    with open(os.path.join(args.out_dir, "train_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("TRAIN_DONE", flush=True)


if __name__ == "__main__":
    main()
