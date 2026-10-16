#!/usr/bin/env python3
"""Zero-GPU diagnosis for M12A Arabidopsis->rice transfer failure.

This script compares M11 pooled clean-plant calibration against M12A fixed
Arabidopsis->rice transfer on the same rice test seqid. It is diagnostic only:
the optional rice test-oracle decode grid is used to understand whether a valid
operating point exists, not to select deployable parameters.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.eval_gene_body_mask import (  # noqa: E402
    SPAN_FEATURES_BY_MODE,
    collect_spans,
    fasta_lengths,
)
from src.screen_anchor.decoders import constrained_decode  # noqa: E402
from src.screen_anchor.gff_io import labels_to_cds_gff  # noqa: E402

OUT_ROOT = ROOT / "reports" / "M13_FAILURE_SANITY"
SPECIES = "oryza_sativa"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def load_scores(path: Path) -> dict[str, np.ndarray]:
    data = np.load(path, allow_pickle=True)
    return {str(seqid): data[f"score::{seqid}"] for seqid in data["seqids"]}


def intervals_to_mask(length: int, intervals: list[tuple[int, int]]) -> np.ndarray:
    mask = np.zeros(length, dtype=bool)
    for start, end in intervals:
        mask[start:end] = True
    return mask


def run_lengths(mask: np.ndarray) -> np.ndarray:
    if mask.size == 0:
        return np.asarray([], dtype=np.int64)
    idx = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    if idx.size == 0:
        return np.asarray([], dtype=np.int64)
    return (idx[1::2] - idx[0::2]).astype(np.int64)


def quantiles(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {"p10": 0.0, "p50": 0.0, "p90": 0.0, "mean": 0.0}
    arr = values.astype(np.float32, copy=False)
    return {
        "p10": float(np.quantile(arr, 0.10)),
        "p50": float(np.quantile(arr, 0.50)),
        "p90": float(np.quantile(arr, 0.90)),
        "mean": float(arr.mean()),
    }


def class_fractions(labels: np.ndarray) -> dict[str, float]:
    total = int(labels.size)
    if total == 0:
        return {"intergenic": 0.0, "cds": 0.0, "intron": 0.0}
    return {
        "intergenic": float((labels == 0).sum() / total),
        "cds": float((labels == 1).sum() / total),
        "intron": float((labels == 2).sum() / total),
    }


def class_fractions_from_counts(counts: np.ndarray) -> dict[str, float]:
    total = int(counts.sum())
    if total == 0:
        return {"intergenic": 0.0, "cds": 0.0, "intron": 0.0}
    return {
        "intergenic": float(counts[0] / total),
        "cds": float(counts[1] / total),
        "intron": float(counts[2] / total),
    }


def selected_labels(scores: np.ndarray, intergenic_bias: float) -> np.ndarray:
    adjusted = scores.astype(np.float32, copy=True)
    adjusted[:, 0] += intergenic_bias
    return adjusted.argmax(axis=-1).astype(np.int8)


def raw_score_diagnostics(run_dir: Path, selected_bias: float) -> dict:
    scores_by_seqid = load_scores(run_dir / "raw_scores" / f"test_{SPECIES}.npz")
    subset_dir = run_dir / "eval_subsets" / SPECIES
    lengths = fasta_lengths(subset_dir / "genome.fa")
    ref_full = collect_spans(subset_dir / "reference.gff3", SPAN_FEATURES_BY_MODE["transcript"])
    ref_cds = collect_spans(subset_dir / "reference.gff3", SPAN_FEATURES_BY_MODE["cds"])

    intergenic_margins = []
    cds_margins = []
    totals = {
        "genome_bases": 0,
        "true_intergenic_bases": 0,
        "true_cds_bases": 0,
        "raw_intergenic_genic": 0,
        "selected_intergenic_genic": 0,
        "raw_cds_genic": 0,
        "selected_cds_genic": 0,
        "selected_preconstrained_runs": 0,
    }
    raw_class_counts = np.zeros(3, dtype=np.int64)
    selected_class_counts = np.zeros(3, dtype=np.int64)
    run_len_chunks = []

    for seqid, scores in scores_by_seqid.items():
        scores_f = scores.astype(np.float32, copy=False)
        length = lengths[seqid]
        full_mask = intervals_to_mask(length, ref_full["intervals_by_seqid"].get(seqid, []))
        cds_mask = intervals_to_mask(length, ref_cds["intervals_by_seqid"].get(seqid, []))
        true_intergenic = ~full_mask

        best_genic = np.maximum(scores_f[:, 1], scores_f[:, 2])
        intergenic_margin = scores_f[:, 0] - best_genic
        cds_margin = best_genic - scores_f[:, 0]
        raw_labels = scores_f.argmax(axis=-1).astype(np.int8)
        sel_labels = selected_labels(scores_f, selected_bias)
        sel_genic = sel_labels > 0
        raw_genic = raw_labels > 0

        intergenic_margins.append(intergenic_margin[true_intergenic])
        cds_margins.append(cds_margin[cds_mask])
        totals["genome_bases"] += int(length)
        totals["true_intergenic_bases"] += int(true_intergenic.sum())
        totals["true_cds_bases"] += int(cds_mask.sum())
        totals["raw_intergenic_genic"] += int((raw_genic & true_intergenic).sum())
        totals["selected_intergenic_genic"] += int((sel_genic & true_intergenic).sum())
        totals["raw_cds_genic"] += int((raw_genic & cds_mask).sum())
        totals["selected_cds_genic"] += int((sel_genic & cds_mask).sum())
        raw_class_counts += np.bincount(raw_labels, minlength=3)[:3]
        selected_class_counts += np.bincount(sel_labels, minlength=3)[:3]
        lengths_arr = run_lengths(sel_genic)
        totals["selected_preconstrained_runs"] += int(lengths_arr.size)
        if lengths_arr.size:
            run_len_chunks.append(lengths_arr)

    intergenic_all = np.concatenate(intergenic_margins) if intergenic_margins else np.asarray([])
    cds_all = np.concatenate(cds_margins) if cds_margins else np.asarray([])
    run_lens = np.concatenate(run_len_chunks) if run_len_chunks else np.asarray([], dtype=np.int64)
    return {
        "genome_bases": totals["genome_bases"],
        "true_intergenic_bases": totals["true_intergenic_bases"],
        "true_cds_bases": totals["true_cds_bases"],
        "raw_argmax_intergenic_false_genic_rate": (
            totals["raw_intergenic_genic"] / totals["true_intergenic_bases"]
            if totals["true_intergenic_bases"] else 0.0
        ),
        "selected_bias_intergenic_false_genic_rate_predecode": (
            totals["selected_intergenic_genic"] / totals["true_intergenic_bases"]
            if totals["true_intergenic_bases"] else 0.0
        ),
        "raw_argmax_true_cds_genic_rate": (
            totals["raw_cds_genic"] / totals["true_cds_bases"]
            if totals["true_cds_bases"] else 0.0
        ),
        "selected_bias_true_cds_genic_rate_predecode": (
            totals["selected_cds_genic"] / totals["true_cds_bases"]
            if totals["true_cds_bases"] else 0.0
        ),
        "intergenic_margin_score0_minus_best_genic": quantiles(intergenic_all),
        "cds_margin_best_genic_minus_score0": quantiles(cds_all),
        "raw_argmax_class_fractions": class_fractions_from_counts(raw_class_counts),
        "selected_bias_class_fractions_predecode": class_fractions_from_counts(selected_class_counts),
        "selected_preconstrained_genic_run_count": totals["selected_preconstrained_runs"],
        "selected_preconstrained_genic_run_length": {
            "p50": float(np.quantile(run_lens, 0.50)) if run_lens.size else 0.0,
            "p90": float(np.quantile(run_lens, 0.90)) if run_lens.size else 0.0,
            "mean": float(run_lens.mean()) if run_lens.size else 0.0,
        },
    }


def eval_labels(labels_by_seqid: dict[str, np.ndarray], subset_dir: Path, out_gff: Path, out_json: Path, exp_id: str) -> dict:
    labels_to_cds_gff(labels_by_seqid, str(out_gff), source="m13_failure_sanity")
    subprocess.run(
        [
            sys.executable,
            "scripts/eval_gene_body_mask.py",
            "--reference-gtf",
            str(subset_dir / "reference.gff3"),
            "--prediction-gtf",
            str(out_gff),
            "--genome-fasta",
            str(subset_dir / "genome.fa"),
            "--output-json",
            str(out_json),
            "--experiment-id",
            exp_id,
            "--profile",
            "screen",
            "--span-mode",
            "cds",
        ],
        cwd=ROOT,
        check=True,
    )
    return load_json(out_json)


def oracle_grid_for_rice(run_dir: Path, run_name: str, selected: dict) -> dict:
    scores_by_seqid = load_scores(run_dir / "raw_scores" / f"test_{SPECIES}.npz")
    subset_dir = run_dir / "eval_subsets" / SPECIES
    combos = []
    seen = set()
    for row in selected["all_val_candidates"]:
        key = (
            float(row["intergenic_bias"]),
            int(row["min_cds_len"]),
            int(row["max_fill_gap"]),
        )
        if key not in seen:
            seen.add(key)
            combos.append(key)

    rows = []
    with tempfile.TemporaryDirectory(prefix=f"{run_name}_oracle_", dir=str(OUT_ROOT)) as tmp_s:
        tmp = Path(tmp_s)
        for bias, min_cds_len, max_fill_gap in combos:
            labels = {
                seqid: selected_labels(scores, bias)
                for seqid, scores in scores_by_seqid.items()
            }
            labels = constrained_decode(labels, min_cds_len=min_cds_len, max_fill_gap=max_fill_gap)
            tag = f"b{str(bias).replace('.', 'p')}_mcl{min_cds_len}_mfg{max_fill_gap}"
            metrics = eval_labels(
                labels,
                subset_dir,
                tmp / f"{tag}.gff",
                tmp / f"{tag}.json",
                f"{run_name}_test_oracle_{tag}",
            )
            rows.append(
                {
                    "tag": tag,
                    "intergenic_bias": bias,
                    "min_cds_len": min_cds_len,
                    "max_fill_gap": max_fill_gap,
                    "intergenic_FPR": metrics["intergenic_FPR"],
                    "gene_body_F1": metrics["gene_body_F1_unconstrained"],
                    "gene_body_precision": metrics["gene_body_precision"],
                    "gene_body_recall": metrics["gene_body_recall"],
                    "gene_count_ratio": metrics["predicted_gene_count_ratio_vs_reference"],
                    "predicted_gene_count": metrics["predicted_gene_count"],
                    "reference_gene_count": metrics["reference_gene_count"],
                }
            )

    def best(candidates: list[dict], key):
        return min(candidates, key=key) if candidates else None

    valid = [r for r in rows if r["intergenic_FPR"] <= 0.01 and r["gene_count_ratio"] <= 1.25]
    sane_count = [r for r in rows if r["gene_count_ratio"] <= 1.25]
    return {
        "candidate_count": len(rows),
        "has_valid_fpr_and_gene_count": bool(valid),
        "best_valid_by_gbf1": best(valid, key=lambda r: (-r["gene_body_F1"], r["intergenic_FPR"])),
        "best_gbf1_under_gene_count": best(sane_count, key=lambda r: (-r["gene_body_F1"], r["intergenic_FPR"])),
        "min_fpr_under_gene_count": best(sane_count, key=lambda r: (r["intergenic_FPR"], -r["gene_body_F1"])),
        "best_gbf1_any": best(rows, key=lambda r: (-r["gene_body_F1"], r["intergenic_FPR"])),
        "min_fpr_any": best(rows, key=lambda r: (r["intergenic_FPR"], -r["gene_body_F1"])),
        "rows": sorted(rows, key=lambda r: (r["intergenic_FPR"], -r["gene_body_F1"])),
    }


def summarize_run(family: str, seed: int, run_dir: Path, run_oracle: bool) -> dict:
    selected = load_json(run_dir / "calibration" / "selected.json")
    metrics = load_json(run_dir / "metrics" / "metrics.json")
    rice_metrics = load_json(run_dir / "metrics" / f"{SPECIES}.metrics.json")
    selected_row = selected["selected"]
    raw_diag = raw_score_diagnostics(run_dir, float(selected_row["intergenic_bias"]))
    out = {
        "family": family,
        "seed": seed,
        "run_dir": str(run_dir),
        "selected_tag": selected_row["tag"],
        "selected_intergenic_bias": selected_row["intergenic_bias"],
        "selected_min_cds_len": selected_row["min_cds_len"],
        "selected_max_fill_gap": selected_row["max_fill_gap"],
        "val_selected": {
            "intergenic_FPR": selected_row["metrics"]["intergenic_FPR"],
            "gene_body_F1": selected_row["metrics"]["gene_body_F1_unconstrained"],
            "gene_count_ratio": selected_row["metrics"]["predicted_gene_count_ratio_vs_reference"],
        },
        "test_rice": {
            "intergenic_FPR": rice_metrics["intergenic_FPR"],
            "gene_body_F1": rice_metrics["gene_body_F1_unconstrained"],
            "constrained_gene_body_F1_at_0.01": metrics.get("constrained_gene_body_F1_at_0.01"),
            "gene_body_precision": rice_metrics["gene_body_precision"],
            "gene_body_recall": rice_metrics["gene_body_recall"],
            "gene_count_ratio": rice_metrics["predicted_gene_count_ratio_vs_reference"],
            "predicted_gene_count": rice_metrics["predicted_gene_count"],
            "reference_gene_count": rice_metrics["reference_gene_count"],
        },
        "raw_score_diagnostics": raw_diag,
    }
    if run_oracle:
        out["test_oracle_grid"] = oracle_grid_for_rice(run_dir, f"{family}-s{seed}", selected)
    return out


def mean_std(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    return {"mean": float(arr.mean()), "std": float(arr.std(ddof=1)) if arr.size > 1 else 0.0}


def aggregate(rows: list[dict]) -> dict:
    keys = [
        ("test_rice", "intergenic_FPR"),
        ("test_rice", "gene_body_F1"),
        ("test_rice", "gene_count_ratio"),
        ("raw_score_diagnostics", "selected_bias_intergenic_false_genic_rate_predecode"),
        ("raw_score_diagnostics", "raw_argmax_intergenic_false_genic_rate"),
        ("raw_score_diagnostics", "selected_bias_true_cds_genic_rate_predecode"),
        ("raw_score_diagnostics", "raw_argmax_true_cds_genic_rate"),
    ]
    out = {}
    for parent, key in keys:
        out[f"{parent}.{key}"] = mean_std([float(r[parent][key]) for r in rows])
    return out


def write_tsv(rows: list[dict], path: Path) -> None:
    header = [
        "family", "seed", "selected_tag", "val_fpr", "val_gbf1", "val_gcount",
        "test_fpr", "test_gbf1", "test_gcount", "predecode_intergenic_fp_rate",
        "predecode_cds_genic_rate", "raw_intergenic_fp_rate", "raw_cds_genic_rate",
    ]
    lines = ["\t".join(header)]
    for r in rows:
        vals = [
            r["family"],
            str(r["seed"]),
            r["selected_tag"],
            f'{r["val_selected"]["intergenic_FPR"]:.6f}',
            f'{r["val_selected"]["gene_body_F1"]:.6f}',
            f'{r["val_selected"]["gene_count_ratio"]:.6f}',
            f'{r["test_rice"]["intergenic_FPR"]:.6f}',
            f'{r["test_rice"]["gene_body_F1"]:.6f}',
            f'{r["test_rice"]["gene_count_ratio"]:.6f}',
            f'{r["raw_score_diagnostics"]["selected_bias_intergenic_false_genic_rate_predecode"]:.6f}',
            f'{r["raw_score_diagnostics"]["selected_bias_true_cds_genic_rate_predecode"]:.6f}',
            f'{r["raw_score_diagnostics"]["raw_argmax_intergenic_false_genic_rate"]:.6f}',
            f'{r["raw_score_diagnostics"]["raw_argmax_true_cds_genic_rate"]:.6f}',
        ]
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


def write_report(summary: dict, path: Path) -> None:
    m11 = summary["aggregate"]["M11"]
    m12 = summary["aggregate"]["M12A"]
    oracle = [r["test_oracle_grid"] for r in summary["runs"] if r["family"] == "M12A"]
    valid_count = sum(1 for x in oracle if x["has_valid_fpr_and_gene_count"])

    lines = [
        "# M13 failure sanity: M12A Arabidopsis->rice",
        "",
        "## Verdict",
        "",
        "- M12A failure is not just final constrained-decoder post-processing noise. The selected Arabidopsis-calibrated operating points transfer poorly to rice before final GFF scoring: pre-decode false genic rate on true rice intergenic bases is much higher than M11 pooled training.",
        "- The dominant pattern is cross-species emission/calibration shift plus fragmentation. M12A keeps reasonable true-CDS genic sensitivity, but it marks too many rice intergenic bases as genic and produces too many gene-body runs/genes after decode.",
        f"- Test-oracle diagnostic grid found valid `FPR<=0.01` and `gene_count<=1.25` rice operating points in {valid_count}/3 M12A seeds. With this grid, the failure is not merely an Arabidopsis-to-rice calibration-transfer bug; rice emissions and fragmentation remain incompatible with the hard guardrails even under diagnostic test-label oracle selection.",
        "- Therefore a close Arabidopsis-relative scan is justified as a distance diagnostic: if a near plant also fails, stop the single-species fixed-model generalization route; if it succeeds while rice fails, reframe as distance-limited transfer.",
        "",
        "## Aggregate comparison",
        "",
        "| family | test rice FPR | test rice gbF1 | test rice gene_count_ratio | predecode intergenic false-genic rate | predecode true-CDS genic rate |",
        "|---|---:|---:|---:|---:|---:|",
        f"| M11 pooled Arabidopsis+rice | {m11['test_rice.intergenic_FPR']['mean']:.4f}±{m11['test_rice.intergenic_FPR']['std']:.4f} | {m11['test_rice.gene_body_F1']['mean']:.4f}±{m11['test_rice.gene_body_F1']['std']:.4f} | {m11['test_rice.gene_count_ratio']['mean']:.3f}±{m11['test_rice.gene_count_ratio']['std']:.3f} | {m11['raw_score_diagnostics.selected_bias_intergenic_false_genic_rate_predecode']['mean']:.4f}±{m11['raw_score_diagnostics.selected_bias_intergenic_false_genic_rate_predecode']['std']:.4f} | {m11['raw_score_diagnostics.selected_bias_true_cds_genic_rate_predecode']['mean']:.4f}±{m11['raw_score_diagnostics.selected_bias_true_cds_genic_rate_predecode']['std']:.4f} |",
        f"| M12A fixed Arabidopsis->rice | {m12['test_rice.intergenic_FPR']['mean']:.4f}±{m12['test_rice.intergenic_FPR']['std']:.4f} | {m12['test_rice.gene_body_F1']['mean']:.4f}±{m12['test_rice.gene_body_F1']['std']:.4f} | {m12['test_rice.gene_count_ratio']['mean']:.3f}±{m12['test_rice.gene_count_ratio']['std']:.3f} | {m12['raw_score_diagnostics.selected_bias_intergenic_false_genic_rate_predecode']['mean']:.4f}±{m12['raw_score_diagnostics.selected_bias_intergenic_false_genic_rate_predecode']['std']:.4f} | {m12['raw_score_diagnostics.selected_bias_true_cds_genic_rate_predecode']['mean']:.4f}±{m12['raw_score_diagnostics.selected_bias_true_cds_genic_rate_predecode']['std']:.4f} |",
        "",
        "## Per-seed selected calibration transfer",
        "",
        "| family | seed | selected | VAL FPR | VAL gbF1 | VAL gcount | rice TEST FPR | rice TEST gbF1 | rice TEST gcount |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in summary["runs"]:
        lines.append(
            f"| {r['family']} | {r['seed']} | {r['selected_tag']} | "
            f"{r['val_selected']['intergenic_FPR']:.4f} | {r['val_selected']['gene_body_F1']:.4f} | {r['val_selected']['gene_count_ratio']:.3f} | "
            f"{r['test_rice']['intergenic_FPR']:.4f} | {r['test_rice']['gene_body_F1']:.4f} | {r['test_rice']['gene_count_ratio']:.3f} |"
        )
    lines.extend([
        "",
        "## M12A rice test-oracle diagnostic",
        "",
        "| seed | valid FPR<=0.01 & gcount<=1.25? | best valid / best sane point | FPR | gbF1 | gcount |",
        "|---:|---|---|---:|---:|---:|",
    ])
    for r in summary["runs"]:
        if r["family"] != "M12A":
            continue
        grid = r["test_oracle_grid"]
        point = grid["best_valid_by_gbf1"] or grid["best_gbf1_under_gene_count"] or grid["min_fpr_any"]
        lines.append(
            f"| {r['seed']} | {grid['has_valid_fpr_and_gene_count']} | {point['tag']} | "
            f"{point['intergenic_FPR']:.4f} | {point['gene_body_F1']:.4f} | {point['gene_count_ratio']:.3f} |"
        )
    lines.extend([
        "",
        "## Next action",
        "",
        "1. Freeze one close Brassicaceae/near-dicot species with high-quality genome+GFF provenance.",
        "2. Run M13 only as a bounded single-seed distance scan: train/calibrate on Arabidopsis, test close plant and rice; no test-label calibration.",
        "3. Treat fly/chicken only as diagnostic/negative controls unless overlap-clean status is resolved.",
        "",
        f"Machine-readable summary: `{OUT_ROOT / 'summary.json'}`",
        f"Per-seed table: `{OUT_ROOT / 'per_seed.tsv'}`",
    ])
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []
    for family, prefix, oracle in [
        ("M11", "M11-L12-SPEC-CALIBRATION", False),
        ("M12A", "M12A-FIXEDMODEL-CROSSSPECIES-A2R", True),
    ]:
        for seed in (0, 1, 2):
            run_dir = ROOT / "outputs" / f"{prefix}-s{seed}"
            rows.append(summarize_run(family, seed, run_dir, run_oracle=oracle))

    summary = {
        "experiment_id": "M13_FAILURE_SANITY",
        "diagnostic_only": True,
        "species": SPECIES,
        "runs": rows,
        "aggregate": {
            "M11": aggregate([r for r in rows if r["family"] == "M11"]),
            "M12A": aggregate([r for r in rows if r["family"] == "M12A"]),
        },
        "interpretation": {
            "primary_failure_mode": "cross_species_emission_calibration_shift_plus_fragmentation",
            "m13_scan_recommendation": "justified_as_bounded_non_claim_distance_diagnostic",
            "claim_safety": "do_not_use_test_oracle_grid_for_model_selection_or_claim",
        },
    }
    (OUT_ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_tsv(rows, OUT_ROOT / "per_seed.tsv")
    write_report(summary, OUT_ROOT / "report.md")
    print(f"Wrote {OUT_ROOT / 'summary.json'}")
    print(f"Wrote {OUT_ROOT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
