"""FP-SEGMENTNT-PROBE-M1 · extract FROZEN SegmentNT base-resolution element logits as input
features for the light probe head. JAX/Haiku (SegmentNT) runs here; the head trains in torch
(separate process) reading the cached .npz — no in-process jax/torch coexistence needed.

Mirrors src/screen_anchor/data.py for split/labels so the probe is a same-budget, same-split
comparison vs the anchor (D.read_fasta / D.build_labels / D.assign_splits — identical
chromosome-level split). Caches per-SEQID full-length (L, F) feature arrays; the head trainer
slices the SAME 2048bp windows from them (window logic stays in the trainer, identical to the
anchor's D.WindowDataset).

TILE extraction (fast + better context): each seqid is processed in TILE_TOKENS-token tiles
(TILE_TOKENS%4==0 for the U-Net's 2 downsample blocks; tile_bp = TILE_TOKENS*6, SegmentNT
trained on 30kb=5000 tokens so a few-kb tile is well within range and gives flanking context).
'N'/non-ACGT -> 'A' before tokenize (SegmentNT 6-mer can't encode 'N'; rare here; noted).

Cache per species: <out-dir>/<species>.npz with
  <seqid> -> (L, F) float16 present-prob per element feature   (stored as feats_<i> + seqid_<i>)
  features (F,) str   config.features order
  Plus split.json-style arrays so the trainer knows train/val/test per seqid.
Run on a GPU node (srun/sbatch), NOT login. WRITE TO SHARED FS (outputs/...), not node-local /tmp.
Usage: python -m src.foundation_probe.extract_segmentnt --species data/m1_screen/<sp> \
         --out-dir outputs/FP-SEGMENTNT-FEATCACHE/<model> [--model segment_nt_multi_species] \
         [--tile-tokens 1000] [--max-seqids 0]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT / "refs/repos/segmentnt-2025"
if __package__ in (None, ""):
    sys.path.insert(0, str(ROOT))
from src.screen_anchor import data as D  # noqa: E402

ACGT = set("ACGT")


def _clean(sub: str) -> str:
    return "".join(c if c in ACGT else "A" for c in sub)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--out-dir", required=True, help="SHARED-FS cache dir; writes <out-dir>/<species>.npz")
    ap.add_argument("--model", default="segment_nt_multi_species")
    ap.add_argument("--tile-tokens", type=int, default=1000, help="tokens per forward tile (%%4==0); bp=tokens*6")
    ap.add_argument("--max-seqids", type=int, default=0, help="0=all; >0 caps for smoke")
    ap.add_argument("--mem-fraction", default="0.9")
    args = ap.parse_args()

    assert args.tile_tokens % 4 == 0, "tile-tokens must be %4==0 (U-Net 2 downsample blocks)"
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = args.mem_fraction
    os.environ.setdefault("TF_GPU_ALLOCATOR", "cuda_malloc_async")
    sys.path.insert(0, str(REPO))

    import jax
    import jax.numpy as jnp
    import haiku as hk
    from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model

    TILE_TOK = args.tile_tokens
    TILE_BP = TILE_TOK * 6
    print(f"[jax] {jax.__version__} devices={jax.devices()}; tile_tokens={TILE_TOK} tile_bp={TILE_BP}", flush=True)
    parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
        model_name=args.model, max_positions=TILE_TOK + 1)
    forward_fn = hk.transform(forward_fn)
    features = list(config.features)
    F = len(features)
    key = jax.random.PRNGKey(0)

    # params passed as ARG (not closure constant) -> avoids the 2.23GB captured-constants blow-up
    apply_jit = jax.jit(forward_fn.apply)

    def infer_tile(seqstr):
        """seqstr length==TILE_BP -> (TILE_BP, F) present-prob float16 (single tile, batch 1)."""
        tok = jnp.asarray([t[1] for t in tokenizer.batch_tokenize([seqstr])], dtype=jnp.int32)
        outs = apply_jit(parameters, key, tok)
        probs = jax.nn.softmax(outs["logits"], axis=-1)[..., -1]   # (1, TILE_BP, F)
        return np.asarray(probs[0], dtype=np.float16)

    sp_path = args.species
    name = os.path.basename(sp_path.rstrip("/"))
    seqs = D.read_fasta(os.path.join(sp_path, "genome.fa"))
    splits = D.assign_splits(list(seqs.keys()))
    seqids = list(seqs.keys())
    if args.max_seqids:
        seqids = seqids[:args.max_seqids]
    print(f"[{name}] seqs={len(seqids)} features={F}", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    save = {"features": np.array(features, dtype=str)}
    seqid_list, split_list = [], []
    for n_done, sid in enumerate(seqids):
        seq = seqs[sid]
        L = len(seq)
        feat = np.zeros((L, F), dtype=np.float16)
        for s in range(0, L, TILE_BP):
            e = min(s + TILE_BP, L)
            sub = _clean(seq[s:e])
            sub = sub + "A" * (TILE_BP - len(sub))         # pad last tile to TILE_BP
            out = infer_tile(sub)                           # (TILE_BP, F)
            feat[s:e] = out[:e - s]
        save[f"feat::{sid}"] = feat
        seqid_list.append(sid); split_list.append(splits[sid])
        if n_done % 5 == 0 or L > 1_000_000:
            print(f"  [{name}] {n_done+1}/{len(seqids)} seqid={sid} L={L}", flush=True)

    save["seqids"] = np.array(seqid_list, dtype=str)
    save["splits"] = np.array(split_list, dtype=str)
    out_path = os.path.join(args.out_dir, f"{name}.npz")
    np.savez_compressed(out_path, **save)
    print(f"[{name}] cached -> {out_path}  seqids={len(seqid_list)} F={F}", flush=True)
    print("EXTRACT_DONE", flush=True)


if __name__ == "__main__":
    main()
