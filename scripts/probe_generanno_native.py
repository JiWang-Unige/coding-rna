#!/usr/bin/env python3
"""GENERANNO-NATIVE-PROBE — run the official GENERanno-1.2b CDS-annotator inference
(NO fine-tuning) on our arabidopsis TEST split, decode per-base CDS predictions with
the OFFICIAL eukaryote pipeline (sliding-window, dual-strand heads, optional official
postprocess), merge the two strand heads into a single gene-body mask, write a CDS GFF
+ test-subset eval_subsets, and call our ruler (eval_gene_body_mask.py --span-mode cds)
to get intergenic_specificity / gene_body_F1_unconstrained on OUR identical ruler.

Decode logic mirrors refs/repos/generanno-2025/src/tasks/downstream/cds_annotation.py:
  - tokenizer(chrs, add_special_tokens=False) -> 1 token == 1 bp (SingleNucleotide k=1)
  - model(input_ids, attention_mask).logits, reshape to (num_heads, valid_len, 2)
  - per-head softmax, overlap-average over the sliding windows, argmax -> per-base CDS(1)/NON(0)
  - head0=positive_strand, head1=negative_strand
  - OFFICIAL eukaryote postprocess (stair-refine + cleanup_short_runs) optional via --postproc
  - merge strands: a base is CDS if EITHER strand head predicts CDS (strand-agnostic ruler)

This is a pure-inference behavior baseline; the screen profile can NEVER claim SOTA.
"""
import argparse
import json
import os
import sys

import numpy as np
import torch
from transformers import AutoTokenizer, AutoModelForTokenClassification

# repo-root on path so we can import our data/gff helpers + the official decode utils
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.screen_anchor import data as D
from src.screen_anchor.gff_io import labels_to_cds_gff
from src.screen_anchor.decoders import constrained_decode

# official postprocess utilities (numba-jit) live in the vendored repo
OFFICIAL_DIR = os.path.join(ROOT, "refs", "repos", "generanno-2025", "src", "tasks", "downstream")
sys.path.insert(0, OFFICIAL_DIR)
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("cds_annotation", os.path.join(OFFICIAL_DIR, "cds_annotation.py"))
_cds = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_cds)  # exposes postprocess_argmax_stair_refine, cleanup_short_binary_runs

LABEL2CHAR = {"CDS": 1, "NON_CODING": 0}


def official_postprocess(argmax_preds, class1_conf,
                         out_shift=128, in_shift=0, stop_run=0, stop_ratio=0.0,
                         min_cds=4, min_gap=4):
    """Apply the OFFICIAL eukaryote postprocess to one head's per-base argmax/conf."""
    refined = _cds.postprocess_argmax_stair_refine(
        class1_confidence=class1_conf, argmax_preds=argmax_preds,
        max_shift=out_shift, inner_shift=in_shift, stop_run=stop_run, stop_ratio=stop_ratio)
    if min_cds > 1 or min_gap > 1:
        refined = _cds.cleanup_short_binary_runs(refined, min_zero_run=min_gap, min_one_run=min_cds)
    return refined


@torch.no_grad()
def annotate_one(seq, model, tokenizer, device, id2label, num_heads, k,
                 ctx_len, overlap_len, micro_bs, dtype,
                 do_postproc, pp_kwargs):
    """Return list[num_heads] of per-base int CDS(1)/NON(0) arrays (len == len(seq)).

    Faithfully mirrors process_sequences_on_gpu chunking + overlap-average decode."""
    pad_id = tokenizer.pad_token_id
    max_char = ctx_len * k
    overlap_char = overlap_len * k
    seq_chunks, seq_masks, seq_chunk_pos = [], [], []
    chr_len = 0

    def add_chunk(chunk_ids):
        pad_len = ctx_len - len(chunk_ids)
        seq_chunks.append(chunk_ids + [pad_id] * pad_len)
        seq_masks.append([1] * len(chunk_ids) + [0] * pad_len)

    seq_work = seq
    if len(seq) < k:  # too short to tokenize
        return [np.zeros(len(seq), dtype=np.int64) for _ in range(num_heads)]
    if len(seq_work) < max_char:
        chrs = seq_work[:len(seq_work) // k * k]
        ids = tokenizer(chrs, add_special_tokens=False)["input_ids"]
        seq_chunk_pos.append((chr_len, 0, len(chrs))); chr_len += len(chrs); add_chunk(ids)
        if len(seq_work) % k != 0:
            chrs = seq_work[len(seq_work) % k:]
            ids = tokenizer(chrs, add_special_tokens=False)["input_ids"]
            seq_chunk_pos.append((chr_len, len(seq_work) % k, len(chrs))); chr_len += len(chrs); add_chunk(ids)
    else:
        while True:
            chrs = seq_work[:max_char]
            ids = tokenizer(chrs, add_special_tokens=False)["input_ids"]
            seq_chunk_pos.append((chr_len, len(seq) - len(seq_work), len(chrs)))
            assert len(chrs) == max_char
            chr_len += len(chrs); add_chunk(ids)
            if len(chrs) == len(seq_work):
                break
            if len(seq_work) - max_char + overlap_char < max_char:
                seq_work = seq_work[-max_char:]
            else:
                seq_work = seq_work[max_char - overlap_char:]

    probs_per_head_chunks = [[] for _ in range(num_heads)]
    total = len(seq_chunks)
    for start in range(0, total, micro_bs):
        end = min(start + micro_bs, total)
        inp = torch.tensor(seq_chunks[start:end], dtype=torch.long, device=device)
        att = torch.tensor(seq_masks[start:end], dtype=torch.long, device=device)
        logits = model(input_ids=inp, attention_mask=att).logits
        probs = logits.softmax(dim=-1).float().cpu()
        for i in range(probs.shape[0]):
            valid_len = int(att[i].sum().item()) * k
            chunk_all = probs[i, :valid_len * num_heads].view(num_heads, valid_len, -1)
            for h in range(num_heads):
                probs_per_head_chunks[h].append(chunk_all[h])

    seq_len = len(seq)
    out_per_head = []
    for h in range(num_heads):
        seq_probs = torch.cat(probs_per_head_chunks[h], dim=0).float()
        final = torch.zeros(seq_len, seq_probs.shape[1], dtype=torch.float32)
        cnt = torch.zeros(seq_len, dtype=torch.long)
        for orig_s, new_s, clen in seq_chunk_pos:
            final[new_s:new_s + clen] += seq_probs[orig_s:orig_s + clen]
            cnt[new_s:new_s + clen] += 1
        final /= cnt.unsqueeze(-1).clamp(min=1)
        argmax = final.argmax(-1).numpy().astype(np.int64)
        if do_postproc and final.shape[1] > 1:
            conf = final[:, 1].numpy()
            argmax = official_postprocess(argmax, conf, **pp_kwargs)
        # map model class id -> {0 NON,1 CDS} via id2label
        out = np.array([LABEL2CHAR[id2label[int(v)]] for v in argmax], dtype=np.int64) \
            if id2label else (argmax != 0).astype(np.int64)
        out_per_head.append(out)
    return out_per_head


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-path", default=os.path.join(ROOT, "data/m1_screen/arabidopsis_thaliana"))
    ap.add_argument("--species-name", default="arabidopsis_thaliana")
    ap.add_argument("--model-name", default="GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--exp-id", default="GENERANNO-NATIVE-PROBE")
    ap.add_argument("--context-length", type=int, default=6000)
    ap.add_argument("--overlap-length", type=int, default=1024)
    ap.add_argument("--micro-bs", type=int, default=1)
    ap.add_argument("--bf16", action="store_true")
    ap.add_argument("--postproc", choices=["on", "off", "constrained"], default="on",
                    help="on = official eukaryote postprocess (true native behavior)")
    ap.add_argument("--limit-seqids", type=int, default=None, help="smoke: only first N test seqids")
    ap.add_argument("--seqids", default=None, help="smoke: explicit comma-separated seqids (override split)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    pred_dir = os.path.join(args.out_dir, "predictions"); os.makedirs(pred_dir, exist_ok=True)
    sub_root = os.path.join(args.out_dir, "eval_subsets"); os.makedirs(sub_root, exist_ok=True)

    seqs = D.read_fasta(os.path.join(args.species_path, "genome.fa"))
    splits = D.assign_splits(list(seqs.keys()))  # i%5==4 -> test (identical ruler as trainer)
    if args.seqids:
        test_ids = [s.strip() for s in args.seqids.split(",") if s.strip() in seqs]
    else:
        test_ids = sorted([sid for sid, t in splits.items() if t == "test"])
        if args.limit_seqids:
            test_ids = test_ids[:args.limit_seqids]
    print(f"[{args.species_name}] total_seqids={len(seqs)} TEST_seqids={len(test_ids)}: {test_ids}", flush=True)
    if not test_ids:
        print("ERROR: no test seqids resolved", flush=True); sys.exit(2)

    dtype = torch.bfloat16 if args.bf16 else torch.float32
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"loading {args.model_name} dtype={dtype} device={device}", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)
    # NOTE: passing dtype=/torch_dtype= as a kwarg makes transformers 4.49 try to JSON-log a
    # config that contains a torch.dtype object -> "Object of type dtype is not JSON serializable".
    # Load in default precision, then cast the module to the target dtype afterwards.
    model = AutoModelForTokenClassification.from_pretrained(
        args.model_name, trust_remote_code=True)
    model.to(device=device, dtype=dtype).eval()
    num_heads = getattr(model, "num_prediction_heads", 2)
    k = getattr(tokenizer, "k", 1) or 1
    id2label = getattr(model.config, "id2label", None)
    print(f"num_heads={num_heads} tokenizer_k={k} id2label={id2label} "
          f"params={sum(p.numel() for p in model.parameters())/1e6:.0f}M", flush=True)

    pp_kwargs = dict(out_shift=128, in_shift=0, stop_run=0, stop_ratio=0.0, min_cds=4, min_gap=4)
    do_pp = args.postproc == "on"

    pred_by_seqid = {}
    for n, sid in enumerate(test_ids, 1):
        seq = seqs[sid].upper()
        heads = annotate_one(seq, model, tokenizer, device, id2label, num_heads, k,
                             args.context_length, args.overlap_length, args.micro_bs, dtype,
                             do_pp, pp_kwargs)
        # merge strands: CDS if EITHER head predicts CDS -> 3class {0 intergenic, 1 CDS}
        merged = np.zeros(len(seq), dtype=np.int8)
        cds_any = np.zeros(len(seq), dtype=bool)
        for h in heads:
            cds_any |= (h != 0)
        merged[cds_any] = D.CLASS_CDS
        pred_by_seqid[sid] = merged
        print(f"  ({n}/{len(test_ids)}) {sid} len={len(seq)} pred_CDS_bases={int(cds_any.sum())} "
              f"({100.0*cds_any.mean():.2f}%)", flush=True)

    if args.postproc == "constrained":
        pred_by_seqid = constrained_decode(pred_by_seqid, min_cds_len=60, max_fill_gap=20)
        print("applied OUR constrained_decode (min_cds_len=60, max_fill_gap=20)", flush=True)
    _src = f"generanno_native_pp-{args.postproc}"
    n_genes = labels_to_cds_gff(pred_by_seqid, os.path.join(pred_dir, f"{args.species_name}.gff"), source=_src)
    print(f"predicted_genes(CDS-runs)={n_genes}", flush=True)

    # write test-subset genome/ref (identical helpers as trainer) for the ruler
    d = os.path.join(sub_root, args.species_name); os.makedirs(d, exist_ok=True)
    D.write_subset_fasta(seqs, test_ids, os.path.join(d, "genome.fa"))
    D.write_subset_gff(os.path.join(args.species_path, "reference.gff3"), test_ids,
                       os.path.join(d, "reference.gff3"))

    summary = {"exp_id": args.exp_id, "species": args.species_name, "model": args.model_name,
               "context_length": args.context_length, "overlap_length": args.overlap_length,
               "postproc": args.postproc, "test_seqids": test_ids, "predicted_genes": n_genes,
               "num_heads": num_heads, "tokenizer_k": k, "bf16": args.bf16}
    with open(os.path.join(args.out_dir, "probe_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("PROBE_DONE", flush=True)


if __name__ == "__main__":
    main()
