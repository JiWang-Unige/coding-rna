#!/usr/bin/env python3
"""Same-panel error analysis for M19 GENERanno and external fixed baselines."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.eval_gene_body_mask import (
    SPAN_FEATURES_BY_MODE,
    collect_spans,
    fasta_lengths,
    intersection_length,
    interval_length,
)


OUT_DIR = ROOT / "reports/M20-SOTA-ERROR-ANALYSIS"
SPECIES = ["arabidopsis_thaliana", "oryza_sativa"]
MODELS = {
    "GENERanno-1.2B-LoRA-s0": {
        "metrics": ROOT / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0/metrics/metrics.json",
        "predictions": ROOT / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0/predictions/{species}.gff",
        "kind": "adapted_pretrained",
    },
    "GENERanno-1.2B-LoRA-s1": {
        "metrics": ROOT / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/metrics/metrics.json",
        "predictions": ROOT / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/predictions/{species}.gff",
        "kind": "adapted_pretrained",
    },
    "Tiberius-angiosperm": {
        "metrics": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/metrics/metrics.json",
        "predictions": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/predictions/{species}.gtf",
        "kind": "released_fixed_model",
    },
    "ANNEVO-Magnoliopsida": {
        "metrics": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/metrics/metrics.json",
        "predictions": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/predictions/{species}.gff",
        "kind": "released_fixed_model",
    },
    "Helixer-land_plant": {
        "metrics": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-HELIXER/metrics/metrics.json",
        "predictions": ROOT / "outputs/M12B-SAMEPANEL-BASELINES-HELIXER/predictions/{species}.gff3",
        "kind": "released_fixed_model",
    },
}
REF_ROOT = ROOT / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/eval_subsets"


def _species_from_key(key: str) -> str | None:
    for species in SPECIES:
        if key.endswith(species):
            return species
    return None


def _load_metrics(model: str, path: Path) -> tuple[dict, list[dict]]:
    data = json.loads(path.read_text())
    aggregate = {
        "model": model,
        "kind": MODELS[model]["kind"],
        "gene_body_F1_unconstrained": data.get("gene_body_F1_unconstrained"),
        "constrained_gene_body_F1": data.get("constrained_gene_body_F1"),
        "gene_body_precision": data.get("gene_body_precision"),
        "gene_body_recall": data.get("gene_body_recall"),
        "intergenic_specificity": data.get("intergenic_specificity"),
        "intergenic_FPR": data.get("intergenic_FPR"),
        "macro_intergenic_specificity": data.get("macro_intergenic_specificity"),
        "predicted_gene_count_ratio_vs_reference": data.get("predicted_gene_count_ratio_vs_reference"),
        "intergenic_guardrail_pass_at_0.01": data.get("intergenic_guardrail_pass_at_0.01"),
    }
    rows = []
    for key, value in data.get("per_species", {}).items():
        species = _species_from_key(key)
        if species is None:
            continue
        rows.append(
            {
                "model": model,
                "kind": MODELS[model]["kind"],
                "species": species,
                "gene_body_F1_unconstrained": value.get("gene_body_F1_unconstrained"),
                "constrained_gene_body_F1": value.get("constrained_gene_body_F1"),
                "gene_body_precision": value.get("gene_body_precision"),
                "gene_body_recall": value.get("gene_body_recall"),
                "intergenic_specificity": value.get("intergenic_specificity"),
                "intergenic_FPR": value.get("intergenic_FPR"),
                "predicted_gene_count": value.get("predicted_gene_count"),
                "reference_gene_count": value.get("reference_gene_count"),
                "predicted_gene_count_ratio_vs_reference": (
                    value.get("predicted_gene_count") / value.get("reference_gene_count")
                    if value.get("predicted_gene_count") is not None and value.get("reference_gene_count")
                    else None
                ),
            }
        )
    return aggregate, rows


def _interval_row(model: str, species: str, pred_path: Path) -> dict:
    ref_path = REF_ROOT / species / "reference.gff3"
    fasta_path = REF_ROOT / species / "genome.fa"
    lengths = fasta_lengths(fasta_path)
    ref = collect_spans(ref_path, SPAN_FEATURES_BY_MODE["cds"])
    pred = collect_spans(pred_path, SPAN_FEATURES_BY_MODE["cds"])
    ref_len = interval_length(ref["intervals_by_seqid"])
    pred_len = interval_length(pred["intervals_by_seqid"])
    overlap = intersection_length(ref["intervals_by_seqid"], pred["intervals_by_seqid"])
    precision = overlap / pred_len if pred_len else 0.0
    recall = overlap / ref_len if ref_len else 0.0
    return {
        "model": model,
        "kind": MODELS[model]["kind"],
        "species": species,
        "genome_bases": sum(lengths.values()),
        "reference_cds_bases": ref_len,
        "predicted_cds_bases": pred_len,
        "overlap_cds_bases": overlap,
        "false_positive_cds_bases": pred_len - overlap,
        "false_negative_cds_bases": ref_len - overlap,
        "cds_base_precision": precision,
        "cds_base_recall": recall,
        "reference_span_groups": ref["span_group_count"],
        "predicted_span_groups": pred["span_group_count"],
        "span_group_ratio_vs_reference": pred["span_group_count"] / ref["span_group_count"]
        if ref["span_group_count"]
        else None,
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _fmt(value: object, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    aggregates = []
    per_species = []
    interval_rows = []
    missing = []
    for model, cfg in MODELS.items():
        if not cfg["metrics"].exists():
            missing.append(str(cfg["metrics"].relative_to(ROOT)))
            continue
        aggregate, species_rows = _load_metrics(model, cfg["metrics"])
        aggregates.append(aggregate)
        per_species.extend(species_rows)
        for species in SPECIES:
            pred_path = Path(str(cfg["predictions"]).format(species=species))
            if not pred_path.exists():
                missing.append(str(pred_path.relative_to(ROOT)))
                continue
            interval_rows.append(_interval_row(model, species, pred_path))

    summary = {
        "exp_id": "M20-SOTA-ERROR-ANALYSIS",
        "species": SPECIES,
        "missing_inputs": missing,
        "aggregate_metrics": aggregates,
        "per_species_metrics": per_species,
        "interval_overlap": interval_rows,
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    _write_csv(OUT_DIR / "aggregate_metrics.csv", aggregates)
    _write_csv(OUT_DIR / "per_species_metrics.csv", per_species)
    _write_csv(OUT_DIR / "interval_overlap.csv", interval_rows)

    ranked = sorted(aggregates, key=lambda row: row.get("constrained_gene_body_F1") or -1.0, reverse=True)
    lines = [
        "# M20-SOTA-ERROR-ANALYSIS",
        "",
        "Same clean-plant panel, same CDS-span evaluator. GENERanno rows are our adapted models; ANNEVO/Tiberius/Helixer rows are released fixed-model baselines.",
        "",
        "## Aggregate Metrics",
        "",
        "| Model | Kind | gbF1 | Precision | Recall | Spec | FPR | Gene count ratio | FPR<=0.01 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in ranked:
        lines.append(
            f"| {row['model']} | {row['kind']} | {_fmt(row['constrained_gene_body_F1'])} | "
            f"{_fmt(row['gene_body_precision'])} | {_fmt(row['gene_body_recall'])} | "
            f"{_fmt(row['intergenic_specificity'])} | {_fmt(row['intergenic_FPR'])} | "
            f"{_fmt(row['predicted_gene_count_ratio_vs_reference'])} | "
            f"{row['intergenic_guardrail_pass_at_0.01']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Tiberius is the strongest released fixed-model comparator under the hard FPR guardrail, but it under-calls gene count relative to reference.",
            "- ANNEVO has the best gbF1 among released fixed baselines on this panel, but aggregate FPR exceeds the `0.01` claim guardrail.",
            "- Helixer strongly over-calls intergenic bases on this panel under the current evaluator, which makes it useful as a practical-specificity contrast.",
            "- GENERanno LoRA is stable across two seeds and keeps FPR under `0.01`, but its remaining weakness is recall/gene recovery rather than specificity. The structured-decoder line should target this exact error mode.",
            "",
            "## Artifacts",
            "",
            "- `summary.json`",
            "- `aggregate_metrics.csv`",
            "- `per_species_metrics.csv`",
            "- `interval_overlap.csv`",
            "",
        ]
    )
    if missing:
        lines.extend(["## Missing Inputs", ""])
        lines.extend(f"- `{item}`" for item in missing)
        lines.append("")
    (OUT_DIR / "report.md").write_text("\n".join(lines))
    print(f"wrote {OUT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
