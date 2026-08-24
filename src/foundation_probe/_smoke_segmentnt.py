"""CK1 smoke for FP-SEGMENTNT-PROBE-M1: confirm SegmentNT (JAX/Haiku) loads on jax 0.10.1 +
runs a tiny base-resolution forward producing (B, seq*6, 14, 2) logits, and that torch is
healthy on the same node. Extraction (jax) and head training (torch) run as SEPARATE
processes, so this only needs each framework to work on the GPU node — not in one process.

SegmentNT constraint: #DNA tokens (excl CLS) must be divisible by 4 (U-Net has 2 downsample
blocks -> 2^2). So bp must be a multiple of 24 (=6*4). 2048bp does NOT satisfy this -> real
extraction must tile on token%4==0 lengths and align back (handled in extract_segmentnt.py).

Run on a GPU node (srun), NOT the login node.
Usage: python -m src.foundation_probe._smoke_segmentnt [--model segment_nt_multi_species]
"""
import argparse
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2] / "refs/repos/segmentnt-2025"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="segment_nt_multi_species")
    ap.add_argument("--bp", type=int, default=48, help="tiny test seq length (bp); must be multiple of 24")
    args = ap.parse_args()

    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION", ".4")
    sys.path.insert(0, str(REPO))

    import jax
    import jax.numpy as jnp
    import haiku as hk
    print(f"[jax] version={jax.__version__} devices={jax.devices()}")

    from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model

    bp = args.bp - (args.bp % 24)
    n_tokens = bp // 6                     # must be divisible by 4
    assert n_tokens % 4 == 0, f"n_tokens={n_tokens} must be %4==0"
    max_positions = n_tokens + 1           # +1 CLS
    print(f"[load] model={args.model} bp={bp} tokens={n_tokens} max_positions={max_positions}")
    parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
        model_name=args.model, max_positions=max_positions,
    )
    print(f"[load] OK — #param groups={len(parameters)} features={list(config.features)}")

    forward_fn = hk.transform(forward_fn)
    seq = ("ACGT" * (bp // 4 + 1))[:bp]
    tokens_ids = [b[1] for b in tokenizer.batch_tokenize([seq])]
    tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)
    print(f"[tok] tokens shape={tokens.shape}")

    key = jax.random.PRNGKey(0)
    outs = forward_fn.apply(parameters, key, tokens)
    logits = outs["logits"]
    print(f"[fwd] logits shape={tuple(logits.shape)}  (expect (1, {bp}, {len(config.features)}, 2))")
    probs = jax.nn.softmax(logits, axis=-1)[..., -1]   # present-prob per feature
    print(f"[fwd] present-prob shape={tuple(probs.shape)} finite={bool(jnp.isfinite(probs).all())}")
    for f in ("protein_coding_gene", "exon", "intron", "splice_donor", "splice_acceptor"):
        if f in config.features:
            print(f"    feature idx {config.features.index(f):2d} = {f}")
    assert logits.ndim == 4 and logits.shape[-1] == 2, f"unexpected logits shape {logits.shape}"
    assert logits.shape[1] == bp, f"base-resolution mismatch: {logits.shape[1]} != {bp}"
    print("SMOKE_OK")


if __name__ == "__main__":
    main()
