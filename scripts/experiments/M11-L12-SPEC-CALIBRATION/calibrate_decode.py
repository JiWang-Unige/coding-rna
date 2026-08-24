#!/usr/bin/env python3
"""Validation-only decode/FPR calibration for M11-L12-SPEC-CALIBRATION.

Consumes raw per-base emissions saved by train_unfreeze_backbone.py, sweeps decode
operating points on VAL only, then applies the selected point once to TEST.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.screen_anchor.decoders import constrained_decode  # noqa: E402
from src.screen_anchor.gff_io import labels_to_cds_gff  # noqa: E402


def parse_csv_floats(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def parse_csv_ints(text: str) -> list[int]:
    return [int(x) for x in text.split(",") if x.strip()]


def tag_for(bias: float, min_cds_len: int, max_fill_gap: int) -> str:
    b = str(bias).replace("-", "m").replace(".", "p")
    return f"b{b}_mcl{min_cds_len}_mfg{max_fill_gap}"


def load_scores(path: Path) -> dict[str, np.ndarray]:
    data = np.load(path, allow_pickle=True)
    seqids = [str(x) for x in data["seqids"]]
    return {sid: data[f"score::{sid}"] for sid in seqids}


def scores_to_labels(scores: dict[str, np.ndarray], intergenic_bias: float) -> dict[str, np.ndarray]:
    out = {}
    for sid, arr in scores.items():
        adjusted = arr.astype(np.float32, copy=True)
        adjusted[:, 0] += intergenic_bias
        out[sid] = adjusted.argmax(axis=-1).astype(np.int8)
    return out


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, cwd=ROOT, check=True)


def evaluate_combo(args, split: str, bias: float, min_cds_len: int, max_fill_gap: int,
                   final_predictions: bool = False) -> dict:
    tag = tag_for(bias, min_cds_len, max_fill_gap)
    if final_predictions:
        pred_dir = Path(args.out_dir) / "predictions"
        metrics_dir = Path(args.out_dir) / "metrics"
        agg_path = metrics_dir / "metrics.json"
        exp_id = args.exp_id
    else:
        pred_dir = Path(args.out_dir) / "calibration" / split / tag / "predictions"
        metrics_dir = Path(args.out_dir) / "calibration" / split / tag / "metrics"
        agg_path = metrics_dir / "metrics.json"
        exp_id = f"{args.exp_id}_{split}_{tag}"
    pred_dir.mkdir(parents=True, exist_ok=True)
    metrics_dir.mkdir(parents=True, exist_ok=True)

    metric_paths = []
    species_for_split = args.val_species if split == "val" else args.test_species
    for sp in species_for_split:
        scores = load_scores(Path(args.out_dir) / "raw_scores" / f"{split}_{sp}.npz")
        pred = scores_to_labels(scores, bias)
        pred = constrained_decode(pred, min_cds_len=min_cds_len, max_fill_gap=max_fill_gap)
        gff_path = pred_dir / f"{sp}.gff"
        n_genes = labels_to_cds_gff(pred, str(gff_path), source=f"m11_cal_{tag}")
        subset_name = sp if split == "test" else f"val_{sp}"
        subset_dir = Path(args.out_dir) / "eval_subsets" / subset_name
        metric_path = metrics_dir / f"{sp}.metrics.json"
        run([
            sys.executable, "scripts/eval_gene_body_mask.py",
            "--reference-gtf", str(subset_dir / "reference.gff3"),
            "--prediction-gtf", str(gff_path),
            "--genome-fasta", str(subset_dir / "genome.fa"),
            "--output-json", str(metric_path),
            "--experiment-id", f"{exp_id}_{sp}",
            "--profile", args.profile,
            "--span-mode", "cds",
        ])
        metric_paths.append(metric_path)
        print(f"[{split} {tag} {sp}] predicted_genes={n_genes}", flush=True)

    run([
        sys.executable, "scripts/aggregate_gene_body_metrics.py",
        "--metrics", *[str(p) for p in metric_paths],
        "--output-json", str(agg_path),
        "--experiment-id", exp_id,
        "--profile", args.profile,
    ])
    return json.loads(agg_path.read_text())


def selection_key(row: dict, target_fpr: float, min_gbf1: float, max_gene_count: float):
    m = row["metrics"]
    fpr = float(m.get("intergenic_FPR", 1.0))
    gbf1 = float(m.get("gene_body_F1_unconstrained", 0.0))
    gcount = float(m.get("predicted_gene_count_ratio_vs_reference", 999.0))
    eligible = fpr <= target_fpr and gbf1 >= min_gbf1 and gcount <= max_gene_count
    if eligible:
        return (0, -gbf1, fpr, abs(1.0 - gcount))
    if gbf1 >= min_gbf1 and gcount <= max_gene_count:
        return (1, fpr, -gbf1, abs(1.0 - gcount))
    return (2, -gbf1, fpr, abs(1.0 - gcount))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-id", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--species", nargs="+", required=True)
    ap.add_argument("--val-species", nargs="*", default=None,
                    help="Optional VAL species list. Defaults to --species for M11 compatibility.")
    ap.add_argument("--test-species", nargs="*", default=None,
                    help="Optional TEST species list. Defaults to --species for M11 compatibility.")
    ap.add_argument("--profile", default="screen")
    ap.add_argument("--intergenic-biases", default="0,0.25,0.5,0.75,1.0,1.25,1.5")
    ap.add_argument("--min-cds-lens", default="60,90")
    ap.add_argument("--max-fill-gaps", default="0,20")
    ap.add_argument("--target-fpr", type=float, default=0.01)
    ap.add_argument("--min-gbf1", type=float, default=0.70)
    ap.add_argument("--max-gene-count-ratio", type=float, default=1.25)
    args = ap.parse_args()
    args.val_species = args.val_species or args.species
    args.test_species = args.test_species or args.species

    out_dir = Path(args.out_dir)
    (out_dir / "calibration").mkdir(parents=True, exist_ok=True)
    rows = []
    for bias in parse_csv_floats(args.intergenic_biases):
        for min_cds_len in parse_csv_ints(args.min_cds_lens):
            for max_fill_gap in parse_csv_ints(args.max_fill_gaps):
                metrics = evaluate_combo(args, "val", bias, min_cds_len, max_fill_gap)
                row = {
                    "tag": tag_for(bias, min_cds_len, max_fill_gap),
                    "intergenic_bias": bias,
                    "min_cds_len": min_cds_len,
                    "max_fill_gap": max_fill_gap,
                    "metrics": metrics,
                }
                rows.append(row)
                print(
                    f"[VAL {row['tag']}] FPR={metrics.get('intergenic_FPR')} "
                    f"gbF1={metrics.get('gene_body_F1_unconstrained')} "
                    f"gcount={metrics.get('predicted_gene_count_ratio_vs_reference')}",
                    flush=True,
                )

    rows.sort(key=lambda r: selection_key(
        r, args.target_fpr, args.min_gbf1, args.max_gene_count_ratio))
    selected = rows[0]
    test_metrics = evaluate_combo(
        args, "test", selected["intergenic_bias"], selected["min_cds_len"],
        selected["max_fill_gap"], final_predictions=True)
    summary = {
        "selection_policy": {
            "target_fpr": args.target_fpr,
            "min_gbf1": args.min_gbf1,
            "max_gene_count_ratio": args.max_gene_count_ratio,
            "criterion": "eligible: max gbF1 then min FPR; fallback: min FPR among sane gbF1/count",
        },
        "selected": selected,
        "test_metrics": test_metrics,
        "all_val_candidates": rows,
    }
    (out_dir / "calibration" / "selected.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"SELECTED {selected['tag']} -> TEST FPR={test_metrics.get('intergenic_FPR')} "
        f"gbF1={test_metrics.get('gene_body_F1_unconstrained')} "
        f"gcount={test_metrics.get('predicted_gene_count_ratio_vs_reference')}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
