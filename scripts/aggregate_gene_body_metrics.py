#!/usr/bin/env python3
"""Aggregate per-species gene-body mask metrics into one lwcr metrics JSON."""

import argparse
import json
from pathlib import Path


def load_metric(path):
    with open(path) as handle:
        return json.load(handle)


def f1(precision, recall):
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def macro_mean(rows, key):
    """Simple unweighted mean over per-species rows (M1-AGGREGATION-GATE-AUDIT).

    Macro weights every species equally, so a single weak species is visible even when
    a base-count-weighted aggregate is dominated by a large genome. Reported alongside
    the base-weighted aggregate and per-species rows; it does NOT drive the gate.
    """
    vals = [r.get(key) for r in rows
            if isinstance(r.get(key), (int, float)) and not isinstance(r.get(key), bool)]
    return sum(vals) / len(vals) if vals else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", nargs="+", required=True, help="Per-species metrics JSON files.")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--profile", default="screen", choices=["smoke", "screen", "full", "scale"])
    parser.add_argument("--smoke-screen-fpr-threshold", type=float, default=0.02)
    parser.add_argument("--full-scale-fpr-threshold", type=float, default=0.01)
    parser.add_argument("--sensitivity-thresholds", default="0.005,0.01,0.02")
    args = parser.parse_args()

    rows = [load_metric(path) for path in args.metrics]
    if not rows:
        raise SystemExit("no metrics provided")

    totals = {
        "genome_bases": 0,
        "reference_gene_body_bases": 0,
        "predicted_gene_body_bases": 0,
        "gene_body_overlap_bases": 0,
        "predicted_intergenic_false_positive_bases": 0,
        "reference_intergenic_bases": 0,
        "reference_gene_count": 0,
        "predicted_gene_count": 0,
        "reference_transcript_count": 0,
        "prediction_transcript_count": 0,
    }
    per_species = {}
    for path, row in zip(args.metrics, rows):
        label = row.get("species") or row.get("data_id") or row.get("experiment_id") or Path(path).stem
        for key in totals:
            totals[key] += row.get(key, 0)
        per_species[label] = {
            "metrics_path": path,
            "intergenic_specificity": row.get("intergenic_specificity"),
            "constrained_gene_body_F1": row.get("constrained_gene_body_F1"),
            "gene_body_F1_unconstrained": row.get("gene_body_F1_unconstrained"),
            "intergenic_FPR": row.get("intergenic_FPR"),
            "gene_body_precision": row.get("gene_body_precision"),
            "gene_body_recall": row.get("gene_body_recall"),
            "reference_gene_count": row.get("reference_gene_count"),
            "predicted_gene_count": row.get("predicted_gene_count"),
        }

    overlap = totals["gene_body_overlap_bases"]
    pred_len = totals["predicted_gene_body_bases"]
    ref_len = totals["reference_gene_body_bases"]
    ref_intergenic = totals["reference_intergenic_bases"]
    pred_only = totals["predicted_intergenic_false_positive_bases"]

    precision = overlap / pred_len if pred_len else 0.0
    recall = overlap / ref_len if ref_len else 0.0
    gene_body_f1 = f1(precision, recall)
    # intergenic = complement of FULL-transcript span (revise-goal 2026-06-11); the per-species
    # evaluator already writes full-transcript-based reference_intergenic_bases + intergenic_FP,
    # so this base-weighted ratio inherits the new ruler. specificity = 1 - FPR (PRIMARY).
    intergenic_fpr = pred_only / ref_intergenic if ref_intergenic else 0.0
    intergenic_specificity = 1.0 - intergenic_fpr
    threshold = args.smoke_screen_fpr_threshold if args.profile in {"smoke", "screen"} else args.full_scale_fpr_threshold
    guardrail_pass = intergenic_fpr <= threshold

    out = {
        "experiment_id": args.experiment_id,
        "profile": args.profile,
        "primary_metric": "intergenic_specificity",
        "semantic_success": True,
        "aggregation_mode": "base-count-weighted_across_species",
        "aggregation_modes_reported": ["base-weighted", "macro", "per-species"],
        "species_count": len(rows),
        "per_species": per_species,
        "intergenic_definition": "complement_of_full_transcript_span_incl_UTR",
        "intergenic_specificity": intergenic_specificity,
        "macro_intergenic_specificity": macro_mean(rows, "intergenic_specificity"),
        "macro_constrained_gene_body_F1": macro_mean(rows, "constrained_gene_body_F1"),
        "macro_gene_body_F1_unconstrained": macro_mean(rows, "gene_body_F1_unconstrained"),
        "macro_intergenic_FPR": macro_mean(rows, "intergenic_FPR"),
        "macro_gene_body_precision": macro_mean(rows, "gene_body_precision"),
        "macro_gene_body_recall": macro_mean(rows, "gene_body_recall"),
        "intergenic_FPR_threshold_used": threshold,
        "intergenic_guardrail_pass": guardrail_pass,
        "constrained_gene_body_F1": gene_body_f1 if guardrail_pass else 0.0,
        "gene_body_F1_unconstrained": gene_body_f1,
        "gene_body_precision": precision,
        "gene_body_recall": recall,
        "intergenic_FPR": intergenic_fpr,
        "predicted_gene_count_ratio_vs_reference": (
            totals["predicted_gene_count"] / totals["reference_gene_count"]
            if totals["reference_gene_count"]
            else 0.0
        ),
        "predicted_transcript_count_ratio_vs_reference": (
            totals["prediction_transcript_count"] / totals["reference_transcript_count"]
            if totals["reference_transcript_count"]
            else 0.0
        ),
        "nucleotide_gene_body_F1_drop_vs_anchor": 0.0,
        **totals,
        "notes": (
            "M1 reproduction aggregate. Counts are base-weighted across species; "
            "nucleotide_gene_body_F1_drop_vs_anchor is set to 0.0 only because no "
            "screen_anchor is frozen yet, and must be recomputed after anchor freeze."
        ),
    }

    for threshold_str in [x.strip() for x in args.sensitivity_thresholds.split(",") if x.strip()]:
        threshold_value = float(threshold_str)
        suffix = threshold_str.rstrip("0").rstrip(".")
        ok = intergenic_fpr <= threshold_value
        out[f"intergenic_guardrail_pass_at_{suffix}"] = ok
        out[f"constrained_gene_body_F1_at_{suffix}"] = gene_body_f1 if ok else 0.0

    output = Path(args.output_json)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        json.dump(out, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
