#!/usr/bin/env python3
"""Build paper-facing same-evaluator comparison tables for M19-era evidence."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]


CLEAN_PLANT_RUNS = [
    {
        "label": "GENERanno 1.2B CDS-preview + our 3-class LoRA, seed0",
        "short": "GENERanno-1.2B-s0",
        "kind": "our adapted/fine-tuned pretrained-CDS backbone",
        "claim_status": "mechanism/challenger; GENERanno provenance overlap unknown",
        "path": "outputs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0/metrics/metrics.json",
    },
    {
        "label": "GENERanno 0.5B base + our 3-class LoRA, seed0",
        "short": "GENERanno-0.5B-base-s0",
        "kind": "our adapted/fine-tuned generic pretrained backbone",
        "claim_status": "mechanism ablation; tests whether generic GENERanno pretraining is sufficient",
        "path": "outputs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0/metrics/metrics.json",
    },
    {
        "label": "GENERanno 1.2B CDS-preview + our 3-class LoRA, M19 seed0",
        "short": "M19-GENERanno-1.2B-s0",
        "kind": "our adapted/fine-tuned pretrained-CDS backbone",
        "claim_status": "pending M19 completion; provenance overlap unknown",
        "path": "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0/metrics/metrics.json",
    },
    {
        "label": "GENERanno 1.2B CDS-preview + our 3-class LoRA, M19 seed1",
        "short": "M19-GENERanno-1.2B-s1",
        "kind": "our adapted/fine-tuned pretrained-CDS backbone",
        "claim_status": "pending M19 completion; provenance overlap unknown",
        "path": "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/metrics/metrics.json",
    },
    {
        "label": "Tiberius angiosperm released model",
        "short": "Tiberius",
        "kind": "released fixed model",
        "claim_status": "same-panel external comparator",
        "path": "outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/metrics/metrics.json",
    },
    {
        "label": "ANNEVO Magnoliopsida released model",
        "short": "ANNEVO",
        "kind": "released fixed model",
        "claim_status": "same-panel external comparator",
        "path": "outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/metrics/metrics.json",
    },
    {
        "label": "Helixer land-plant released model",
        "short": "Helixer",
        "kind": "released fixed model",
        "claim_status": "same-panel external comparator",
        "path": "outputs/M12B-SAMEPANEL-BASELINES-HELIXER/metrics/metrics.json",
    },
]


BROAD_PANEL_RUNS = [
    {
        "label": "ANNEVO released clade-matched models",
        "short": "ANNEVO-M17",
        "kind": "released fixed/clade-matched models",
        "claim_status": "diagnostic; overlap/provenance caveats vary by species",
        "path": "outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-ANNEVO/metrics/metrics.json",
    },
    {
        "label": "Tiberius released model",
        "short": "Tiberius-M17",
        "kind": "released fixed model",
        "claim_status": "diagnostic; strong specificity but under-calls genes",
        "path": "outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-TIBERIUS/metrics/metrics.json",
    },
    {
        "label": "Helixer released lineage models",
        "short": "Helixer-M17",
        "kind": "released fixed/clade-matched models",
        "claim_status": "diagnostic; high FPR on plant/gallus rows",
        "path": "outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-HELIXER/metrics/metrics.json",
    },
]


FIELDS = [
    "short",
    "label",
    "kind",
    "claim_status",
    "status",
    "gene_body_F1_unconstrained",
    "constrained_gene_body_F1_at_0.005",
    "constrained_gene_body_F1_at_0.01",
    "constrained_gene_body_F1_at_0.02",
    "intergenic_specificity",
    "intergenic_FPR",
    "macro_intergenic_specificity",
    "predicted_gene_count",
    "reference_gene_count",
    "predicted_gene_count_ratio_vs_reference",
    "utility_note",
    "metrics_path",
]


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def fmt(value: Any, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def utility_note(metrics: dict[str, Any] | None) -> str:
    if metrics is None:
        return "pending"
    fpr = float(metrics.get("intergenic_FPR", 1.0))
    gcount = float(metrics.get("predicted_gene_count_ratio_vs_reference", 999.0))
    notes = []
    if fpr <= 0.01:
        notes.append("FPR<=0.01")
    elif fpr <= 0.02:
        notes.append("FPR 0.01-0.02")
    else:
        notes.append("FPR>0.02")
    if gcount < 0.75:
        notes.append("under-calls genes")
    elif gcount > 1.25:
        notes.append("over-calls genes")
    else:
        notes.append("gene count sane")
    return "; ".join(notes)


def gene_count_ratio(metrics: dict[str, Any]) -> float | None:
    ratio = metrics.get("predicted_gene_count_ratio_vs_reference")
    if ratio is not None:
        return float(ratio)
    pred = metrics.get("predicted_gene_count")
    ref = metrics.get("reference_gene_count")
    if pred is None or ref in (None, 0):
        return None
    return float(pred) / float(ref)


def constrained_at_001(metrics: dict[str, Any]) -> float | None:
    if metrics.get("constrained_gene_body_F1_at_0.01") is not None:
        return metrics.get("constrained_gene_body_F1_at_0.01")
    fpr = metrics.get("intergenic_FPR")
    gbf1 = metrics.get("gene_body_F1_unconstrained")
    if fpr is None or gbf1 is None:
        return None
    return float(gbf1) if float(fpr) <= 0.01 else 0.0


def row_from_spec(spec: dict[str, str], root: Path) -> dict[str, Any]:
    metrics_path = root / spec["path"]
    metrics = load_json(metrics_path)
    row = {k: None for k in FIELDS}
    row.update({
        "short": spec["short"],
        "label": spec["label"],
        "kind": spec["kind"],
        "claim_status": spec["claim_status"],
        "status": "available" if metrics is not None else "pending",
        "metrics_path": str(metrics_path.relative_to(root)),
        "utility_note": utility_note(metrics),
    })
    if metrics is None:
        return row
    for key in FIELDS:
        if key in metrics:
            row[key] = metrics[key]
    row["predicted_gene_count_ratio_vs_reference"] = gene_count_ratio(metrics)
    return row


def species_rows(spec: dict[str, str], root: Path) -> list[dict[str, Any]]:
    metrics = load_json(root / spec["path"])
    if not metrics:
        return []
    out = []
    for key, value in sorted(metrics.get("per_species", {}).items()):
        species = key.split("_", 1)[1] if "_" in key else key
        ratio = gene_count_ratio(value)
        value_for_note = dict(value)
        value_for_note["predicted_gene_count_ratio_vs_reference"] = ratio
        row = {
            "model": spec["short"],
            "species": species,
            "gene_body_F1_unconstrained": value.get("gene_body_F1_unconstrained"),
            "constrained_gene_body_F1_at_0.01": constrained_at_001(value),
            "intergenic_specificity": value.get("intergenic_specificity"),
            "intergenic_FPR": value.get("intergenic_FPR"),
            "predicted_gene_count": value.get("predicted_gene_count"),
            "reference_gene_count": value.get("reference_gene_count"),
            "predicted_gene_count_ratio_vs_reference": ratio,
            "utility_note": utility_note(value_for_note),
        }
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def md_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> list[str]:
    lines = [
        "| " + " | ".join(title for title, _ in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        cells = []
        for _, key in columns:
            digits = 3 if "count_ratio" in key else 4
            cells.append(fmt(row.get(key), digits=digits))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def main() -> int:
    out_dir = ROOT / "reports" / "M19-COMPARABILITY-EVIDENCE"
    out_dir.mkdir(parents=True, exist_ok=True)

    clean_rows = [row_from_spec(spec, ROOT) for spec in CLEAN_PLANT_RUNS]
    broad_rows = [row_from_spec(spec, ROOT) for spec in BROAD_PANEL_RUNS]
    clean_species = [r for spec in CLEAN_PLANT_RUNS for r in species_rows(spec, ROOT)]
    broad_species = [r for spec in BROAD_PANEL_RUNS for r in species_rows(spec, ROOT)]

    write_csv(out_dir / "clean_plant_aggregate.csv", clean_rows, FIELDS)
    write_csv(out_dir / "broad_panel_aggregate.csv", broad_rows, FIELDS)
    species_fields = [
        "model",
        "species",
        "gene_body_F1_unconstrained",
        "constrained_gene_body_F1_at_0.01",
        "intergenic_specificity",
        "intergenic_FPR",
        "predicted_gene_count",
        "reference_gene_count",
        "predicted_gene_count_ratio_vs_reference",
        "utility_note",
    ]
    write_csv(out_dir / "clean_plant_per_species.csv", clean_species, species_fields)
    write_csv(out_dir / "broad_panel_per_species.csv", broad_species, species_fields)

    payload = {
        "clean_plant_aggregate": clean_rows,
        "clean_plant_per_species": clean_species,
        "broad_panel_aggregate": broad_rows,
        "broad_panel_per_species": broad_species,
        "interpretation": {
            "main_message": (
                "GENERanno 1.2B is specificity-controlled and gene-count sane on the "
                "clean plant panel, but released Tiberius/ANNEVO/Helixer still set a "
                "higher clean-plant gbF1 frontier. Its current paper role is "
                "pretrained-CDS adaptation/challenger evidence unless provenance clears."
            ),
            "claim_boundary": (
                "Because public GENERanno sources do not expose a complete species/accession "
                "training manifest, Arabidopsis/rice results are not clean held-out claim "
                "evidence for GENERanno-based models."
            ),
        },
    }
    (out_dir / "comparison_tables.json").write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# M19 Comparability Evidence",
        "",
        "Scope: same-evaluator, same clean-plant panel comparison for paper-facing utility tables. "
        "Rows marked pending will be filled automatically after M19 seed metrics exist.",
        "",
        "## Clean Plant Aggregate",
    ]
    lines.extend(md_table(clean_rows, [
        ("Model", "short"),
        ("Type", "kind"),
        ("Status", "status"),
        ("gbF1", "gene_body_F1_unconstrained"),
        ("gbF1@0.005", "constrained_gene_body_F1_at_0.005"),
        ("gbF1@0.01", "constrained_gene_body_F1_at_0.01"),
        ("gbF1@0.02", "constrained_gene_body_F1_at_0.02"),
        ("Spec", "intergenic_specificity"),
        ("FPR", "intergenic_FPR"),
        ("Macro spec", "macro_intergenic_specificity"),
        ("Gene count ratio", "predicted_gene_count_ratio_vs_reference"),
        ("Utility", "utility_note"),
    ]))
    lines.extend([
        "",
        "## Clean Plant Per Species",
    ])
    lines.extend(md_table(clean_species, [
        ("Model", "model"),
        ("Species", "species"),
        ("gbF1", "gene_body_F1_unconstrained"),
        ("gbF1@0.01", "constrained_gene_body_F1_at_0.01"),
        ("Spec", "intergenic_specificity"),
        ("FPR", "intergenic_FPR"),
        ("Gene count ratio", "predicted_gene_count_ratio_vs_reference"),
        ("Utility", "utility_note"),
    ]))
    lines.extend([
        "",
        "## Broad Diagnostic Aggregate (M17)",
    ])
    lines.extend(md_table(broad_rows, [
        ("Model", "short"),
        ("Type", "kind"),
        ("gbF1", "gene_body_F1_unconstrained"),
        ("gbF1@0.01", "constrained_gene_body_F1_at_0.01"),
        ("Spec", "intergenic_specificity"),
        ("FPR", "intergenic_FPR"),
        ("Macro spec", "macro_intergenic_specificity"),
        ("Gene count ratio", "predicted_gene_count_ratio_vs_reference"),
        ("Utility", "utility_note"),
    ]))
    lines.extend([
        "",
        "## Interpretation",
        "",
        "- M18 GENERanno 1.2B is not a random baseline: it is FPR-valid and gene-count sane on clean plants, unlike the 0.5B base result.",
        "- Clean-plant released callers still define the high-gbF1 frontier: ANNEVO/Tiberius/Helixer are around 0.922-0.927 gbF1; Tiberius is the closest practical comparator because it also passes FPR<=0.01, but it under-calls genes.",
        "- Current GENERanno evidence should be written as pretrained-CDS backbone adaptation/comparability, not clean no-overlap held-out SOTA, until provenance clears or a cleaner species panel is selected.",
        "",
        "## Artifacts",
        "",
        "- `clean_plant_aggregate.csv`",
        "- `clean_plant_per_species.csv`",
        "- `broad_panel_aggregate.csv`",
        "- `broad_panel_per_species.csv`",
        "- `comparison_tables.json`",
    ])
    (out_dir / "comparison_tables.md").write_text("\n".join(lines) + "\n")
    print(out_dir / "comparison_tables.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
