"""M10-GENERANNO-LORA-3C: LoRA-tune GENERanno encoder with our intron-aware 3-class head.

GENERanno's released CDS annotator is a binary CDS/non-coding token classifier.  Its native
head has no intron concept, so intron-sized gaps split one real gene into many predicted genes.
This trainer keeps the pretrained GENERanno encoder, discards the binary CDS head, adds LoRA
adapters to attention projections, and trains the existing per-base 3-class FP-aware head:

  0 intergenic, 1 CDS, 2 gene-body non-CDS (mostly intron/UTR under the current ruler).

Screen-only.  It writes CDS GFF predictions plus eval_subsets in the same layout as the M9
NT-v2 unfreeze trainer.
"""

import argparse
import json
import os
import re
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if __package__ in (None, ""):
    sys.path.insert(0, ROOT)

from src.foundation_probe.train_probe_head import _ConvLSTMHead, fp_penalty, macro_f1  # noqa: E402
from src.screen_anchor import data as D  # noqa: E402
from src.screen_anchor.gff_io import labels_to_cds_gff  # noqa: E402


MODEL = "GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
_NONACGT = re.compile(r"[^ACGT]")


def gene_body_tversky_loss(emissions, Y, alpha=0.75, beta=0.25, ignore=-100, eps=1e-6):
    """Binary gene-body Tversky loss over 3-class emissions.

    CE still handles the 3-class CDS/non-CDS/intergenic labels. This auxiliary term targets
    the publication guardrail directly: gene-body probability on true intergenic bases is FP,
    and intergenic probability on true gene-body bases is FN.
    """
    import torch

    mask = Y != ignore
    if not mask.any():
        return emissions.sum() * 0.0
    p = torch.softmax(emissions, dim=-1)
    p_gene = 1.0 - p[..., 0]
    y_gene = (Y > 0).to(dtype=p_gene.dtype)
    valid = mask.to(dtype=p_gene.dtype)
    tp = (p_gene * y_gene * valid).sum()
    fp = (p_gene * (1.0 - y_gene) * valid).sum()
    fn = ((1.0 - p_gene) * y_gene * valid).sum()
    score = (tp + eps) / (tp + alpha * fp + beta * fn + eps)
    return 1.0 - score


def _clean(seq):
    return _NONACGT.sub("A", seq.upper())


class WindowTokDataset:
    """Uniform full windows only.  GENERanno uses k-mer tokenization and predicts k bases/token."""

    def __init__(
        self,
        seqs,
        labels,
        seqids,
        tokenizer,
        window,
        k,
        sample_fraction=1.0,
        seed=0,
        max_windows=None,
        pretokenize=False,
    ):
        import torch  # noqa: F401

        self.seqs = seqs
        self.labels = labels
        self.tok = tokenizer
        self.window = window
        self.k = k
        self.index = []
        for sid in seqids:
            for start in range(0, len(seqs[sid]), window):
                end = start + window
                if end <= len(seqs[sid]):
                    self.index.append((sid, start, end))
        if sample_fraction < 1.0 and self.index:
            rng = np.random.default_rng(seed)
            n = max(1, int(round(len(self.index) * sample_fraction)))
            chosen = rng.choice(len(self.index), size=n, replace=False)
            self.index = [self.index[i] for i in sorted(chosen)]
        if max_windows is not None and max_windows > 0:
            self.index = self.index[:max_windows]
        self._id_cache = None
        if pretokenize:
            self._id_cache = []
            for sid, start, end in self.index:
                seq = _clean(self.seqs[sid][start:end])
                self._id_cache.append(_tokenize_window(self.tok, seq, self.window, self.k).cpu())

    def __len__(self):
        return len(self.index)

    def __getitem__(self, idx):
        import torch

        sid, start, end = self.index[idx]
        if self._id_cache is None:
            seq = _clean(self.seqs[sid][start:end])
            ids = _tokenize_window(self.tok, seq, self.window, self.k)
        else:
            ids = self._id_cache[idx]
        y = self.labels[sid][start:end].astype(np.int64)
        return ids, torch.ones_like(ids), torch.from_numpy(y)


def _tokenize_window(tokenizer, seq, window, k):
    ids = tokenizer(seq, add_special_tokens=False, return_tensors="pt")["input_ids"][0]
    if len(ids) * k != window:
        raise ValueError(f"token/base alignment failed: {len(ids)} tokens x k={k} != {window}")
    return ids


def collate(batch):
    import torch

    ids = torch.stack([row[0] for row in batch])
    attn = torch.stack([row[1] for row in batch])
    labels = torch.stack([row[2] for row in batch])
    return ids, attn, labels


def _pooled(species, split):
    seqs, labels, ids = {}, {}, []
    for name, sp in species.items():
        for sid, tag in sp["splits"].items():
            if tag != split:
                continue
            key = f"{name}::{sid}"
            seqs[key] = sp["seqs"][sid]
            labels[key] = sp["labels"][sid]
            ids.append(key)
    return seqs, labels, ids


def _write_subsets(sp, ids, root, name):
    out = os.path.join(root, name)
    os.makedirs(out, exist_ok=True)
    D.write_subset_fasta(sp["seqs"], ids, os.path.join(out, "genome.fa"))
    D.write_subset_gff(os.path.join(sp["sp_path"], "reference.gff3"), ids, os.path.join(out, "reference.gff3"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--exp-id", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model-name", default=MODEL)
    parser.add_argument("--window", type=int, default=6144, help="bp window; must be divisible by model k")
    parser.add_argument("--sample-fraction", type=float, default=0.1)
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--patience", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--predict-batch", type=int, default=1)
    parser.add_argument("--lr", type=float, default=8e-4)
    parser.add_argument("--lora-lr", type=float, default=2e-5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--class-weighting", default="sqrt_inv", choices=["none", "inv", "sqrt_inv"])
    parser.add_argument("--loss", default="fp_aware", choices=["ce", "fp_aware", "gb_tversky"])
    parser.add_argument("--fp-lambda", type=float, default=1.0)
    parser.add_argument("--gb-tversky-lambda", type=float, default=1.0)
    parser.add_argument("--gb-tversky-alpha", type=float, default=0.75)
    parser.add_argument("--gb-tversky-beta", type=float, default=0.25)
    parser.add_argument(
        "--decoder",
        default="none",
        choices=["none", "crf"],
        help="Optional structured decoder trained on top of per-base emissions. crf uses LinearChainCRFVec.",
    )
    parser.add_argument(
        "--crf-aux-ce",
        type=float,
        default=1.0,
        help="Auxiliary class-weighted CE weight when --decoder crf is enabled.",
    )
    parser.add_argument("--postproc", default="constrained", choices=["none", "constrained"])
    parser.add_argument("--min-cds-len", type=int, default=60)
    parser.add_argument("--max-fill-gap", type=int, default=20)
    parser.add_argument("--lora-r", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument("--lora-dropout", type=float, default=0.05)
    parser.add_argument("--lora-targets", default="q_proj,k_proj,v_proj,o_proj")
    parser.add_argument("--model-task", default="auto", choices=["auto", "token_classification", "masked_lm"])
    parser.add_argument("--attn-implementation", default="sdpa", choices=["eager", "sdpa", "flash_attention_2"])
    parser.add_argument("--bf16", action="store_true")
    parser.add_argument("--limit-train-windows", type=int, default=None)
    parser.add_argument("--limit-val-windows", type=int, default=None)
    parser.add_argument(
        "--pretokenize-windows",
        action="store_true",
        help="Precompute token ids for train/val windows before training to avoid per-batch tokenizer overhead.",
    )
    parser.add_argument(
        "--train-progress-every",
        type=int,
        default=0,
        help="Print batch-level training progress every N batches; 0 disables.",
    )
    parser.add_argument(
        "--eval-batch",
        type=int,
        default=None,
        help="Validation batch size. Defaults to --batch-size; can be larger because validation is no-grad.",
    )
    parser.add_argument(
        "--eval-progress-every",
        type=int,
        default=0,
        help="Print batch-level validation progress every N batches; 0 disables.",
    )
    parser.add_argument("--limit-test-seqids", type=int, default=None)
    parser.add_argument(
        "--limit-predict-windows",
        type=int,
        default=None,
        help="Optional smoke/debug cap on prediction windows per predict_scores call.",
    )
    parser.add_argument(
        "--save-raw-scores",
        action="store_true",
        help="Save pre-postproc per-base logits for VAL and TEST under raw_scores/*.npz for no-leak calibration.",
    )
    args = parser.parse_args()

    import torch
    from peft import LoraConfig, get_peft_model
    from torch.utils.data import DataLoader
    from transformers import AutoConfig, AutoModelForMaskedLM, AutoModelForTokenClassification, AutoTokenizer

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.bfloat16 if args.bf16 else torch.float32
    os.makedirs(args.out_dir, exist_ok=True)
    pred_dir = os.path.join(args.out_dir, "predictions")
    subset_dir = os.path.join(args.out_dir, "eval_subsets")
    os.makedirs(pred_dir, exist_ok=True)
    os.makedirs(subset_dir, exist_ok=True)

    print(f"loading {args.model_name} attn={args.attn_implementation} dtype={dtype} device={device}", flush=True)
    tok = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)
    cfg = AutoConfig.from_pretrained(args.model_name, trust_remote_code=True)
    arch = ",".join(getattr(cfg, "architectures", []) or [])
    model_task = args.model_task
    if model_task == "auto":
        model_task = "masked_lm" if "MaskedLM" in arch else "token_classification"
    loader = AutoModelForMaskedLM if model_task == "masked_lm" else AutoModelForTokenClassification
    full_model = loader.from_pretrained(
        args.model_name,
        trust_remote_code=True,
        attn_implementation=args.attn_implementation,
    )
    k = getattr(full_model, "k", getattr(full_model.config, "k", None))
    if k is None:
        probe_ids = tok("A" * args.window, add_special_tokens=False, return_tensors="pt")["input_ids"][0]
        if len(probe_ids) <= 0 or args.window % len(probe_ids) != 0:
            raise ValueError(
                f"cannot infer token/base ratio for {args.model_name}: "
                f"{args.window} bp -> {len(probe_ids)} tokens"
            )
        k = args.window // len(probe_ids)
    k = int(k)
    if args.window % k != 0:
        raise ValueError(f"--window {args.window} must be divisible by k={k}")

    backbone = full_model.model
    del full_model
    try:
        backbone.gradient_checkpointing_enable()
        print("gradient_checkpointing enabled", flush=True)
    except (AttributeError, ValueError) as exc:
        print(f"gradient_checkpointing unavailable: {exc}", flush=True)
    backbone.to(device=device, dtype=dtype)
    for param in backbone.parameters():
        param.requires_grad_(False)

    lora_cfg = LoraConfig(
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        target_modules=[x.strip() for x in args.lora_targets.split(",") if x.strip()],
        lora_dropout=args.lora_dropout,
        bias="none",
    )
    backbone = get_peft_model(backbone, lora_cfg)
    if hasattr(backbone, "enable_input_require_grads"):
        backbone.enable_input_require_grads()
    hidden_dim = int(backbone.config.hidden_size)
    trainable_lora = sum(p.numel() for p in backbone.parameters() if p.requires_grad)
    print(
        f"backbone=GENERanno task={model_task} arch={arch or 'unknown'} hidden={hidden_dim} "
        f"layers={backbone.config.num_hidden_layers} k={k} trainable_lora_params={trainable_lora}",
        flush=True,
    )

    probe = tok("A" * args.window, add_special_tokens=False)
    if len(probe["input_ids"]) * k != args.window:
        raise ValueError(f"token alignment check failed: {len(probe['input_ids'])} tokens x{k} != {args.window}")
    print(f"token alignment: {len(probe['input_ids'])} tokens x{k} = {args.window} bp", flush=True)
    with torch.no_grad():
        probe_ids = _tokenize_window(tok, "A" * args.window, args.window, k).unsqueeze(0).to(device)
        probe_attn = torch.ones_like(probe_ids)
        probe_out = backbone(input_ids=probe_ids, attention_mask=probe_attn)
        if not hasattr(probe_out, "last_hidden_state"):
            raise ValueError(f"{args.model_name} backbone output has no last_hidden_state")
        hidden_tokens = int(probe_out.last_hidden_state.shape[1])
        if hidden_tokens * k != args.window:
            raise ValueError(f"hidden/token alignment failed: {hidden_tokens} tokens x k={k} != {args.window}")
        print(f"backbone probe OK: hidden_shape={tuple(probe_out.last_hidden_state.shape)}", flush=True)

    species = {}
    for sp_path in args.species:
        name = os.path.basename(sp_path.rstrip("/"))
        seqs = D.read_fasta(os.path.join(sp_path, "genome.fa"))
        lengths = {sid: len(seq) for sid, seq in seqs.items()}
        labels = D.build_labels(os.path.join(sp_path, "reference.gff3"), lengths)
        splits = D.assign_splits(list(seqs.keys()))
        species[name] = {"sp_path": sp_path, "seqs": seqs, "labels": labels, "splits": splits}
        print(
            f"[{name}] seqs={len(seqs)} splits="
            f"{ {split: sum(1 for val in splits.values() if val == split) for split in ('train', 'val', 'test')} }",
            flush=True,
        )

    tr_s, tr_l, tr_ids = _pooled(species, "train")
    va_s, va_l, va_ids = _pooled(species, "val")
    train_ds = WindowTokDataset(
        tr_s,
        tr_l,
        tr_ids,
        tok,
        args.window,
        k,
        sample_fraction=args.sample_fraction,
        seed=args.seed,
        max_windows=args.limit_train_windows,
        pretokenize=args.pretokenize_windows,
    )
    val_ds = WindowTokDataset(
        va_s,
        va_l,
        va_ids,
        tok,
        args.window,
        k,
        max_windows=args.limit_val_windows,
        pretokenize=args.pretokenize_windows,
    )
    print(
        f"train_windows={len(train_ds)} val_windows={len(val_ds)} "
        f"pretokenize={args.pretokenize_windows}",
        flush=True,
    )
    train_dl = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, collate_fn=collate)
    eval_batch = args.batch_size if args.eval_batch is None else max(1, args.eval_batch)
    val_dl = DataLoader(val_ds, batch_size=eval_batch, shuffle=False, collate_fn=collate)

    cls = np.zeros(D.NUM_CLASSES, dtype=np.int64)
    for sid in tr_ids:
        cls += np.bincount(tr_l[sid], minlength=D.NUM_CLASSES)
    if args.class_weighting == "none" or (cls == 0).any():
        weight = None
    else:
        raw = (1.0 / cls) if args.class_weighting == "inv" else (1.0 / np.sqrt(cls))
        weight = torch.tensor(raw / raw.mean(), dtype=torch.float32, device=device)
    print(f"class_counts={cls.tolist()} class_weighting={args.class_weighting}", flush=True)

    head = _ConvLSTMHead(hidden_dim, n_classes=D.NUM_CLASSES).to(device)
    decoder = None
    if args.decoder == "crf":
        from src.screen_anchor.decoders import LinearChainCRFVec

        decoder = LinearChainCRFVec(D.NUM_CLASSES).to(device)
    opt = torch.optim.Adam(
        [
            {"params": head.parameters(), "lr": args.lr},
            *([{"params": decoder.parameters(), "lr": args.lr}] if decoder is not None else []),
            {"params": [p for p in backbone.parameters() if p.requires_grad], "lr": args.lora_lr},
        ]
    )
    ce = torch.nn.CrossEntropyLoss(ignore_index=-100, weight=weight)
    use_amp = device == "cuda" and args.bf16

    def forward(ids, attn):
        out = backbone(input_ids=ids, attention_mask=attn)
        hidden = out.last_hidden_state
        per_base = hidden.repeat_interleave(k, dim=1)
        if per_base.shape[1] != args.window:
            per_base = per_base[:, : args.window]
        return head(per_base.transpose(1, 2))

    def compute_loss(logits, y):
        aux = ce(logits.reshape(-1, D.NUM_CLASSES), y.reshape(-1))
        if args.loss == "fp_aware":
            aux = aux + args.fp_lambda * fp_penalty(logits, y)
        elif args.loss == "gb_tversky":
            aux = aux + args.gb_tversky_lambda * gene_body_tversky_loss(
                logits,
                y,
                alpha=args.gb_tversky_alpha,
                beta=args.gb_tversky_beta,
            )
        if decoder is not None:
            mask = y != -100
            return decoder.nll(logits, y, mask) + args.crf_aux_ce * aux
        return aux

    def evaluate():
        backbone.eval()
        head.eval()
        if decoder is not None:
            decoder.eval()
        conf = np.zeros((D.NUM_CLASSES, D.NUM_CLASSES), dtype=np.int64)
        val_start = time.time()
        with torch.no_grad():
            for batch_idx, (ids, attn, y) in enumerate(val_dl, start=1):
                ids = ids.to(device)
                attn = attn.to(device)
                y_dev = y.to(device)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                    logits = forward(ids, attn)
                    if decoder is None:
                        pred = logits.argmax(-1).cpu().numpy()
                    else:
                        pred = decoder.viterbi(logits, y_dev != -100).cpu().numpy()
                yy = y.numpy()
                for true, got in zip(yy.reshape(-1), pred.reshape(-1)):
                    if true != -100:
                        conf[int(true), int(got)] += 1
                if args.eval_progress_every > 0 and (
                    batch_idx == 1 or batch_idx % args.eval_progress_every == 0 or batch_idx == len(val_dl)
                ):
                    elapsed = time.time() - val_start
                    print(
                        f"eval_progress batch={batch_idx}/{len(val_dl)} elapsed_s={elapsed:.1f}",
                        flush=True,
                    )
        return macro_f1(conf)

    best_f1, best_head, best_lora, best_decoder, bad = -1.0, None, None, None, 0
    for epoch in range(1, args.epochs + 1):
        backbone.train()
        head.train()
        if decoder is not None:
            decoder.train()
        total, batches = 0.0, 0
        epoch_start = time.time()
        for ids, attn, y in train_dl:
            ids = ids.to(device)
            attn = attn.to(device)
            y = y.to(device)
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                loss = compute_loss(forward(ids, attn), y)
            loss.backward()
            opt.step()
            total += float(loss.item())
            batches += 1
            if args.train_progress_every > 0 and (
                batches == 1 or batches % args.train_progress_every == 0 or batches == len(train_dl)
            ):
                elapsed = time.time() - epoch_start
                print(
                    f"train_progress epoch={epoch} batch={batches}/{len(train_dl)} "
                    f"loss={float(loss.item()):.4f} elapsed_s={elapsed:.1f}",
                    flush=True,
                )
        val_f1, per_class = evaluate()
        print(
            f"epoch {epoch}: train_loss={total / max(batches, 1):.4f} "
            f"val_macroF1={val_f1:.4f} per_class={[round(x, 4) for x in per_class]}",
            flush=True,
        )
        if val_f1 > best_f1:
            best_f1, bad = val_f1, 0
            best_head = {k_: v.detach().cpu().clone() for k_, v in head.state_dict().items()}
            best_lora = {k_: v.detach().cpu().clone() for k_, v in backbone.state_dict().items() if "lora_" in k_}
            best_decoder = (
                None
                if decoder is None
                else {k_: v.detach().cpu().clone() for k_, v in decoder.state_dict().items()}
            )
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop epoch={epoch} best_val_macroF1={best_f1:.4f}", flush=True)
                break
    if best_head is not None:
        head.load_state_dict(best_head)
    if best_lora is not None:
        current = backbone.state_dict()
        current.update(best_lora)
        backbone.load_state_dict(current)
    if decoder is not None and best_decoder is not None:
        decoder.load_state_dict(best_decoder)

    backbone.eval()
    head.eval()
    if decoder is not None:
        decoder.eval()
    source = f"generanno_lora3c_r{args.lora_r}_{args.loss}_{args.decoder}_pp-{args.postproc}"
    raw_score_dir = os.path.join(args.out_dir, "raw_scores")
    if args.save_raw_scores:
        os.makedirs(raw_score_dir, exist_ok=True)

    def predict_scores(sp, seqids):
        scores = {sid: np.zeros((len(sp["seqs"][sid]), D.NUM_CLASSES), dtype=np.float16) for sid in seqids}
        windows = [
            (sid, start, start + args.window)
            for sid in seqids
            for start in range(0, len(sp["seqs"][sid]), args.window)
            if start + args.window <= len(sp["seqs"][sid])
        ]
        if args.limit_predict_windows is not None and args.limit_predict_windows > 0:
            windows = windows[: args.limit_predict_windows]
        print(f"predict_scores seqids={len(seqids)} windows={len(windows)} batch={args.predict_batch}", flush=True)
        with torch.no_grad():
            for offset in range(0, len(windows), args.predict_batch):
                chunk = windows[offset : offset + args.predict_batch]
                ids = []
                for sid, start, end in chunk:
                    ids.append(_tokenize_window(tok, _clean(sp["seqs"][sid][start:end]), args.window, k))
                ids_t = torch.stack(ids).to(device)
                attn_t = torch.ones_like(ids_t)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                    logits = forward(ids_t, attn_t).float().cpu().numpy()
                for i, (sid, start, end) in enumerate(chunk):
                    scores[sid][start:end, :] = logits[i].astype(np.float16)
                step = offset // max(args.predict_batch, 1) + 1
                if offset == 0 or step % 500 == 0 or offset + args.predict_batch >= len(windows):
                    done = min(offset + args.predict_batch, len(windows))
                    print(f"predict_scores progress {done}/{len(windows)}", flush=True)
        return scores

    def labels_from_scores(score_by_seqid, intergenic_bias=0.0):
        pred = {}
        for sid, score in score_by_seqid.items():
            adjusted = score.astype(np.float32, copy=True)
            adjusted[:, 0] += intergenic_bias
            if decoder is None:
                pred[sid] = adjusted.argmax(axis=-1).astype(np.int8)
                continue
            labels = np.zeros(adjusted.shape[0], dtype=np.int8)
            decode_bases = adjusted.shape[0]
            if args.limit_predict_windows is not None and args.limit_predict_windows > 0:
                decode_bases = min(decode_bases, args.limit_predict_windows * args.window)
            starts = [start for start in range(0, decode_bases, args.window) if start + args.window <= decode_bases]
            for offset in range(0, len(starts), max(args.predict_batch, 1)):
                chunk = starts[offset : offset + max(args.predict_batch, 1)]
                emissions_np = np.stack([adjusted[start : start + args.window] for start in chunk], axis=0)
                emissions = torch.from_numpy(emissions_np).to(device=device, dtype=torch.float32)
                mask = torch.ones((len(chunk), args.window), dtype=torch.bool, device=device)
                paths = decoder.viterbi(emissions, mask).cpu().numpy().astype(np.int8)
                for i, start in enumerate(chunk):
                    labels[start : start + args.window] = paths[i]
                step = offset // max(args.predict_batch, 1) + 1
                if offset == 0 or step % 100 == 0 or offset + max(args.predict_batch, 1) >= len(starts):
                    done = min(offset + max(args.predict_batch, 1), len(starts))
                    print(f"crf_decode progress {sid} {done}/{len(starts)} windows", flush=True)
            pred[sid] = labels
        return pred

    def save_scores(split, name, score_by_seqid):
        if not args.save_raw_scores:
            return
        payload = {"seqids": np.array(sorted(score_by_seqid), dtype=str)}
        payload.update({f"score::{sid}": score_by_seqid[sid] for sid in sorted(score_by_seqid)})
        np.savez(os.path.join(raw_score_dir, f"{split}_{name}.npz"), **payload)

    summary = {
        "exp_id": args.exp_id,
        "model": args.model_name,
        "model_task": model_task,
        "seed": args.seed,
        "best_val_macro_f1": best_f1,
        "window": args.window,
        "k": k,
        "sample_fraction": args.sample_fraction,
        "lora_r": args.lora_r,
        "lora_alpha": args.lora_alpha,
        "lora_targets": args.lora_targets,
        "loss": args.loss,
        "fp_lambda": args.fp_lambda,
        "gb_tversky_lambda": args.gb_tversky_lambda,
        "gb_tversky_alpha": args.gb_tversky_alpha,
        "gb_tversky_beta": args.gb_tversky_beta,
        "decoder": args.decoder,
        "crf_aux_ce": args.crf_aux_ce,
        "bf16": args.bf16,
        "eval_batch": eval_batch,
        "device": device,
        "raw_scores_saved": bool(args.save_raw_scores),
    }
    for name, sp in species.items():
        if args.save_raw_scores:
            val_ids = [sid for sid, tag in sp["splits"].items() if tag == "val"]
            if val_ids:
                val_scores = predict_scores(sp, val_ids)
                save_scores("val", name, val_scores)
                _write_subsets(sp, val_ids, subset_dir, f"val_{name}")
        test_ids = [sid for sid, tag in sp["splits"].items() if tag == "test"]
        if args.limit_test_seqids is not None and args.limit_test_seqids > 0:
            test_ids = test_ids[: args.limit_test_seqids]
        if not test_ids:
            continue
        test_scores = predict_scores(sp, test_ids)
        save_scores("test", name, test_scores)
        out = labels_from_scores(test_scores)
        if args.postproc == "constrained":
            from src.screen_anchor.decoders import constrained_decode

            out = constrained_decode(out, min_cds_len=args.min_cds_len, max_fill_gap=args.max_fill_gap)
        genes = labels_to_cds_gff(out, os.path.join(pred_dir, f"{name}.gff"), source=source)
        _write_subsets(sp, test_ids, subset_dir, name)
        summary.setdefault("species", {})[name] = {"test_seqids": test_ids, "predicted_genes": genes}
        print(f"[{name}] test_seqids={len(test_ids)} predicted_genes={genes}", flush=True)

    with open(os.path.join(args.out_dir, "train_summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print("TRAIN_DONE", flush=True)


if __name__ == "__main__":
    main()
