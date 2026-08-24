#!/usr/bin/env python3
"""Test-label oracle decode sweep for M18 raw scores.

This is a NON-CLAIM diagnostic: it intentionally tunes decode parameters on
test labels to determine whether a failed species is a calibration problem or
an emission/coherence problem.
"""
from __future__ import annotations

import argparse
import json
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


def metric_row(metrics: dict, tag: str, bias: float, min_cds_len: int, max_fill_gap: int) -> dict:
    return {
        "tag": tag,
        "intergenic_bias": bias,
        "min_cds_len": min_cds_len,
        "max_fill_gap": max_fill_gap,
        "intergenic_FPR": metrics.get("intergenic_FPR"),
        "intergenic_specificity": metrics.get("intergenic_specificity"),
        "gene_body_F1_unconstrained": metrics.get("gene_body_F1_unconstrained"),
        "constrained_gene_body_F1": metrics.get("constrained_gene_body_F1"),
        "predicted_gene_count": metrics.get("predicted_gene_count"),
        "reference_gene_count": metrics.get("reference_gene_count"),
        "predicted_gene_count_ratio_vs_reference": metrics.get("predicted_gene_count_ratio_vs_reference"),
    }


def best_rows(rows: list[dict]) -> dict:
    def f(x, default):
        return default if x is None else float(x)

    valid_001 = [
        r for r in rows
        if f(r["intergenic_FPR"], 1.0) <= 0.01
        and f(r["predicted_gene_count_ratio_vs_reference"], 999.0) <= 1.25
    ]
    valid_002 = [
        r for r in rows
        if f(r["intergenic_FPR"], 1.0) <= 0.02
        and f(r["predicted_gene_count_ratio_vs_reference"], 999.0) <= 1.25
    ]
    sane_count = [r for r in rows if f(r["predicted_gene_count_ratio_vs_reference"], 999.0) <= 1.25]
    return {
        "best_valid_fpr_0p01": max(valid_001, key=lambda r: f(r["gene_body_F1_unconstrained"], -1.0), default=None),
        "best_valid_fpr_0p02": max(valid_002, key=lambda r: f(r["gene_body_F1_unconstrained"], -1.0), default=None),
        "lowest_fpr_sane_count": min(sane_count, key=lambda r: f(r["intergenic_FPR"], 1.0), default=None),
        "best_gbf1_sane_count": max(sane_count, key=lambda r: f(r["gene_body_F1_unconstrained"], -1.0), default=None),
        "best_gbf1_overall": max(rows, key=lambda r: f(r["gene_body_F1_unconstrained"], -1.0), default=None),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-run-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--species", nargs="+", default=["gallus_gallus"])
    ap.add_argument("--intergenic-biases", default="0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8")
    ap.add_argument("--min-cds-lens", default="60,90,120")
    ap.add_argument("--max-fill-gaps", default="0,20")
    args = ap.parse_args()

    source = Path(args.source_run_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_results = {
        "source_run_dir": str(source),
        "claim_status": "NON-CLAIM test-label oracle diagnostic",
        "species": {},
    }

    for sp in args.species:
        scores = load_scores(source / "raw_scores" / f"test_{sp}.npz")
        species_dir = out_dir / sp
        pred_dir = species_dir / "predictions"
        metrics_dir = species_dir / "metrics"
        pred_dir.mkdir(parents=True, exist_ok=True)
        metrics_dir.mkdir(parents=True, exist_ok=True)
        subset_dir = source / "eval_subsets" / sp
        rows = []
        for bias in parse_csv_floats(args.intergenic_biases):
            for min_cds_len in parse_csv_ints(args.min_cds_lens):
                for max_fill_gap in parse_csv_ints(args.max_fill_gaps):
                    tag = tag_for(bias, min_cds_len, max_fill_gap)
                    pred = scores_to_labels(scores, bias)
                    pred = constrained_decode(pred, min_cds_len=min_cds_len, max_fill_gap=max_fill_gap)
                    gff_path = pred_dir / f"{tag}.gff"
                    labels_to_cds_gff(pred, str(gff_path), source=f"m18_oracle_{tag}")
                    metric_path = metrics_dir / f"{tag}.metrics.json"
                    run([
                        sys.executable, "scripts/eval_gene_body_mask.py",
                        "--reference-gtf", str(subset_dir / "reference.gff3"),
                        "--prediction-gtf", str(gff_path),
                        "--genome-fasta", str(subset_dir / "genome.fa"),
                        "--output-json", str(metric_path),
                        "--experiment-id", f"M18-ORACLE-{sp}_{tag}",
                        "--profile", "screen",
                        "--span-mode", "cds",
                    ])
                    metrics = json.loads(metric_path.read_text())
                    row = metric_row(metrics, tag, bias, min_cds_len, max_fill_gap)
                    rows.append(row)
                    print(
                        f"[{sp} {tag}] FPR={row['intergenic_FPR']} "
                        f"gbF1={row['gene_body_F1_unconstrained']} "
                        f"gcount={row['predicted_gene_count_ratio_vs_reference']}",
                        flush=True,
                    )
        rows_sorted = sorted(rows, key=lambda r: (
            999.0 if r["intergenic_FPR"] is None else float(r["intergenic_FPR"]),
            -(0.0 if r["gene_body_F1_unconstrained"] is None else float(r["gene_body_F1_unconstrained"])),
        ))
        all_results["species"][sp] = {
            "summary": best_rows(rows),
            "rows": rows_sorted,
        }

    (out_dir / "oracle_per_species_calibration.json").write_text(json.dumps(all_results, indent=2) + "\n")
    lines = [
        "# M18 Per-Species Oracle Calibration Diagnostic",
        "",
        "NON-CLAIM: this uses test labels to choose decode parameters. It is only a failure-mode diagnostic.",
        "",
    ]
    for sp, block in all_results["species"].items():
        lines.append(f"## {sp}")
        for name, row in block["summary"].items():
            lines.append(f"- {name}: `{row['tag']}` FPR={row['intergenic_FPR']:.4f} "
                         f"gbF1={row['gene_body_F1_unconstrained']:.4f} "
                         f"gcount={row['predicted_gene_count_ratio_vs_reference']:.3f}" if row else f"- {name}: none")
        lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
