"""Emit chromosome-level split manifests for the leakage gate (scripts/check_data.py).
One JSONL record per seqid: {"id": "<species>::<seqid>", "split": "train|val|test"}.
Because the split is assigned per whole seqid, train vs held-out IDs must be disjoint —
check_data verifies this deterministically before any training.

Usage: python scripts/screen_anchor_make_manifest.py --species data/m1_screen/A data/m1_screen/B --out-dir outputs/<exp>/splits
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.screen_anchor import data as D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    train, heldout, rows = [], [], []
    for sp_path in args.species:
        name = os.path.basename(sp_path.rstrip("/"))
        seqs = D.read_fasta(os.path.join(sp_path, "genome.fa"))
        splits = D.assign_splits(list(seqs.keys()))
        for sid, tag in splits.items():
            rec = {"id": f"{name}::{sid}", "split": tag}
            rows.append(rec)
            (train if tag == "train" else heldout).append(rec)

    with open(os.path.join(args.out_dir, "train.jsonl"), "w") as fh:
        for r in train:
            fh.write(json.dumps(r) + "\n")
    with open(os.path.join(args.out_dir, "heldout.jsonl"), "w") as fh:
        for r in heldout:
            fh.write(json.dumps(r) + "\n")
    with open(os.path.join(args.out_dir, "all_splits.jsonl"), "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    n = {k: sum(1 for r in rows if r["split"] == k) for k in ("train", "val", "test")}
    print(f"manifest: {len(rows)} seqids -> {n}  (train={len(train)}, heldout={len(heldout)})")


if __name__ == "__main__":
    main()
