#!/usr/bin/env python3
"""Generate a small set of result visualizations from docs/06_results_log.md and reports/.

Outputs are written to the repository-level `figure/` directory.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path("/home/users/j/jwang/coding-rna")
RESULTS_MD = BASE_DIR / "docs/06_results_log.md"
REPORTS_DIR = BASE_DIR / "reports"
FIG_DIR = BASE_DIR / "figure"
NT_EXPERIMENT_HINTS = (
    "fp-segnt",
    "fp-fragfix",
    "ta-foundation-decoder-m4",
    "ta-coherence-fix-m5",
    "ta-fragfix-sweep-m6",
    "reanchor-heldout-m7",
    "tb-gbf1-multiclass-m8",
    "tb-unfreeze-backbone-m9",
    "tb-unfreeze-backbone-m9-deep",
    "m9-deep",
    "m10-m9l12",
    "generanno-native-probe",
    "generanno-lora",
    "m10-generanno-lora-3c",
    "nt-v2",
)
NT_ROW_HINTS = (
    "fp-segnt",
    "fp-fragfix",
    "segnt",
    "fragfix",
    "constr",
    "unfreeze",
    "frozen",
    "candidate",
    "anchor",
    "mc-candidate",
    "3c-candidate",
    "raw-dna",
    "generanno",
    "native",
    "lora",
    "l0",
    "l1",
    "l2",
    "l4",
    "l6",
    "l8",
    "l10",
    "l12",
    "s0",
    "s1",
    "s2",
    "s3",
    "s4",
)


def parse_markdown_row(line: str) -> List[str]:
    """Split a markdown table row into stripped cells."""
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def parse_float(value: str) -> float:
    """Parse a numeric-like metric value from markdown/table text."""
    if value is None:
        return float("nan")
    text = value.strip().replace("`", "")
    if not text or text in {"-", "N/A", "na", "NA", "NaN", "nan"}:
        return float("nan")
    # Keep the first number-like token; this also handles accidental punctuation.
    m = re.search(r"[-+]?(\d*\.\d+|\d+)([eE][-+]?\d+)?", text)
    if not m:
        return float("nan")
    try:
        return float(m.group(0))
    except ValueError:
        return float("nan")


def parse_std_from_value(value: str) -> float:
    """Parse the second component from a value like ``mean ± std``."""
    if value is None:
        return float("nan")
    text = str(value).strip().replace("`", "")
    if not text:
        return float("nan")
    if "±" in text:
        m = re.search(r"±\s*([-+]?\d*(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                return float("nan")
    m = re.search(r"\+\s*/\s*-\s*([-+]?\d*(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return float("nan")
    return float("nan")


def parse_markdown_table(lines: Sequence[str], start: int) -> Tuple[Optional[pd.DataFrame], int]:
    """Parse a markdown table starting at `start`.

    Returns (DataFrame, next_index).
    """
    idx = start
    while idx < len(lines) and not lines[idx].strip():
        idx += 1

    if idx >= len(lines) or not lines[idx].lstrip().startswith("|"):
        return None, idx

    header = parse_markdown_row(lines[idx])
    if idx + 1 >= len(lines):
        return None, idx + 1

    sep = lines[idx + 1].strip()
    parts = [p.strip() for p in sep.strip("|").split("|")]
    if len(parts) < 2 or not all(re.fullmatch(r":?-{3,}:?", p) for p in parts if p):
        return None, idx + 1

    rows: List[List[str]] = []
    row_idx = idx + 2
    while row_idx < len(lines):
        if not lines[row_idx].lstrip().startswith("|"):
            break
        row_cells = parse_markdown_row(lines[row_idx])
        if len(row_cells) != len(header):
            # Skip malformed rows and stop only if extremely inconsistent.
            row_idx += 1
            continue
        rows.append(row_cells)
        row_idx += 1

    if not rows:
        return None, row_idx

    df = pd.DataFrame(rows, columns=header)
    return df, row_idx


def is_nt_experiment(name: str) -> bool:
    text = normalize_key(name)
    if any(token in text for token in NT_EXPERIMENT_HINTS):
        return True
    return bool(re.search(r"\b(tb-unfreeze-backbone-m9|m10-m9l12|generanno|nt-v2|fp-segnt|fp-fragfix)\b", text))


def is_nt_row_entry(experiment: str, row: str) -> bool:
    exp = normalize_key(experiment)
    text = normalize_key(row)
    if not text:
        return False
    if any(token in text for token in NT_ROW_HINTS):
        return True
    if "l" in text and re.search(r"\bl[0-9]+\b", text):
        return True
    if "s" in text and re.search(r"\bs[0-9]+\b", text):
        return True
    return is_nt_experiment(exp) and "candidate" in text


def infer_nt_family(experiment: str, row: str) -> str:
    text = normalize_key(f"{experiment} {row}")
    if "fp-segnt" in text:
        return "NT-SEGNT"
    if "fp-fragfix" in text:
        return "NT-FRAGFIX"
    if "tb-gbf1-multiclass-m8" in text or "multiclass" in text:
        return "M8-MULTICLASS"
    if "m10-m9l12" in text:
        return "M10-M9L12"
    if "tb-unfreeze-backbone-m9" in text or "m9-deep" in text or "unfreeze" in text:
        return "M9-UNFREEZE"
    if "m7" in text and "reanchor" in text:
        return "M7-HELDOUT"
    if "m7" in text:
        return "M7-CLEAN"
    if "generanno" in text:
        return "GENERANNO"
    return "NT-OTHER"


def parse_nt_candidate_summary_rows(experiment: str, table_df: pd.DataFrame) -> List[Dict[str, float]]:
    first_col = table_df.columns[0]
    normalized_cols = [normalize_key(c) for c in table_df.columns]
    candidate_cols = [col for col, norm in zip(table_df.columns, normalized_cols) if "candidate" in norm]
    if not candidate_cols:
        return []
    candidate_col = candidate_cols[0]

    spec = float("nan")
    f1 = float("nan")
    fpr = float("nan")

    for _, row in table_df.iterrows():
        metric = normalize_key(row[first_col])
        value = row.get(candidate_col, "")
        if "intergenic_specificity" in metric:
            spec = parse_float(str(value))
        elif "macro_intergenic_specificity" in metric:
            continue
        elif "gene_body_f1" in metric:
            f1 = parse_float(str(value))
        elif "intergenic_fpr" in metric:
            fpr = parse_float(str(value))

    if np.isnan(spec) and np.isnan(f1) and np.isnan(fpr):
        return []

    return [
        {
            "experiment": experiment,
            "run": "candidate-summary",
            "f1": f1,
            "specificity": spec,
            "fpr": fpr,
            "f1_std": float("nan"),
            "spec_std": float("nan"),
            "fpr_std": float("nan"),
            "family": infer_nt_family(experiment, "candidate-summary"),
        }
    ]


def parse_nt_rows_from_markdown_table(experiment: str, table_df: pd.DataFrame) -> List[Dict[str, float]]:
    first_col = table_df.columns[0]
    f1_col = choose_metric_column(table_df.columns, "f1")
    spec_col = choose_metric_column(table_df.columns, "specificity")
    fpr_col = choose_metric_column(table_df.columns, "fpr")

    rows: List[Dict[str, float]] = []
    for _, row in table_df.iterrows():
        run = str(row[first_col]).strip()
        if not is_nt_row_entry(experiment, run):
            continue

        f1 = parse_float(str(row.get(f1_col, ""))) if f1_col else float("nan")
        spec = parse_float(str(row.get(spec_col, ""))) if spec_col else float("nan")
        fpr = parse_float(str(row.get(fpr_col, ""))) if fpr_col else float("nan")
        if np.isnan(spec) and not np.isnan(fpr):
            spec = max(0.0, 1.0 - fpr)
        if np.isnan(f1) and np.isnan(spec):
            continue

        rows.append(
            {
                "experiment": experiment,
                "run": run,
                "f1": f1,
                "specificity": spec,
                "fpr": fpr,
                "f1_std": parse_std_from_value(str(row.get(f1_col, ""))) if f1_col else float("nan"),
                "spec_std": parse_std_from_value(str(row.get(spec_col, ""))) if spec_col else float("nan"),
                "fpr_std": parse_std_from_value(str(row.get(fpr_col, ""))) if fpr_col else float("nan"),
                "family": infer_nt_family(experiment, run),
            }
        )

    if rows:
        return rows
    return parse_nt_candidate_summary_rows(experiment, table_df)


def parse_results_log_nt_points(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    records: List[Dict[str, float]] = []
    current = None

    idx = 0
    while idx < len(lines):
        line = lines[idx]
        if line.startswith("## Result:"):
            current = line.split(":", 1)[1].strip()
            idx += 1
            continue

        if not current:
            idx += 1
            continue

        if not is_nt_experiment(current):
            idx += 1
            continue

        if not line.lstrip().startswith("|"):
            idx += 1
            continue

        table_df, next_idx = parse_markdown_table(lines, idx)
        if table_df is None:
            idx += 1
            continue

        for rec in parse_nt_rows_from_markdown_table(current, table_df):
            records.append(rec)

        idx = next_idx

    return pd.DataFrame.from_records(records)


def normalize_key(name: str) -> str:
    return re.sub(r"\s+", " ", str(name).strip().lower())


def choose_metric_column(columns: Sequence[str], mode: str) -> Optional[str]:
    """Choose the most relevant metric column among candidates."""
    cols = [c for c in columns if c and isinstance(c, str)]
    lc = [normalize_key(c) for c in cols]
    if mode == "f1":
        candidates = []
        for name, norm in zip(cols, lc):
            if "gbf1" not in norm and "gene_body_f1" not in norm:
                continue
            score = 0
            if "constrained" in norm:
                score += 10
            m = re.search(r"@([0-9.]+)", norm)
            if m:
                # prefer tighter threshold first
                score += 20 - int(float(m.group(1)) * 100) if float(m.group(1)) <= 0.1 else 5
            if norm == "gbf1":
                score += 1
            candidates.append((score, name))
        if candidates:
            return sorted(candidates, key=lambda x: x[0], reverse=True)[0][1]
        return None
    if mode == "specificity":
        # Prefer direct intergenic specificity
        direct = [name for name in cols if "specificity" in normalize_key(name)]
        if direct:
            direct = sorted(
                direct,
                key=lambda c: 0 if normalize_key(c) == "specificity" else 1,
            )
            return direct[0]
        direct = [name for name in cols if "specificity" in name.lower()]
        return direct[0] if direct else None
    if mode == "fpr":
        direct = [name for name in cols if re.search(r"\bfpr\b", normalize_key(name))]
        return direct[0] if direct else None
    return None


def pick_aggregate_row(df: pd.DataFrame) -> Optional[Dict[str, str]]:
    if df.empty or " " not in df.columns:
        # keep generic fallback to first column
        first_col = df.columns[0]
    else:
        first_col = df.columns[0]

    candidates = []
    for _, row in df.iterrows():
        key = str(row[first_col]).strip().lower()
        if "aggregate" in key:
            candidates.append(row)

    if not candidates:
        return None

    # Prefer rows that clearly look like TEST/Split aggregate over verbose descriptions.
    candidates = [
        r
        for r in candidates
        if re.match(r"\s*(aggregate|test aggregate|split / species|all|clean|broad)", str(r[first_col]).strip().lower())
    ] or candidates
    return dict(candidates[0])


def parse_results_log_aggregates(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    records = []
    current = None

    for idx, line in enumerate(lines):
        if line.startswith("## Result:"):
            current = line.split(":", 1)[1].strip()
            continue

        if not current:
            continue

        if line.startswith("### Metrics"):
            table_df, next_idx = parse_markdown_table(lines, idx + 1)
            if table_df is None:
                continue
            row = pick_aggregate_row(table_df)
            if row is None:
                continue

            f1_col = choose_metric_column(table_df.columns, "f1")
            spec_col = choose_metric_column(table_df.columns, "specificity")
            fpr_col = choose_metric_column(table_df.columns, "fpr")

            f1 = parse_float(str(row.get(f1_col, ""))) if f1_col else float("nan")
            spec = parse_float(str(row.get(spec_col, ""))) if spec_col else float("nan")
            fpr = parse_float(str(row.get(fpr_col, ""))) if fpr_col else float("nan")
            if np.isnan(spec) and not np.isnan(fpr):
                spec = max(0.0, 1.0 - fpr)

            if np.isnan(f1) and np.isnan(spec):
                continue

            first_col_key = table_df.columns[0]
            records.append(
                {
                    "experiment": current,
                    "f1": f1,
                    "specificity": spec,
                    "fpr": fpr,
                    "aggregate_row": str(row.get(first_col_key, "")),
                }
            )
    return pd.DataFrame.from_records(records)


def plot_clean_plant_frontier() -> None:
    agg_path = (
        REPORTS_DIR / "M19-COMPARABILITY-EVIDENCE" / "clean_plant_aggregate.csv"
    )
    if not agg_path.exists():
        print(f"[warn] missing {agg_path}, skip frontier plot")
        return

    df = pd.read_csv(agg_path)
    df["kind"] = df["kind"].fillna("unknown")
    df["model_short"] = df["short"].astype(str)
    df["f1"] = df["constrained_gene_body_F1_at_0.01"].fillna(df["gene_body_F1_unconstrained"])
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["f1", "intergenic_specificity"])

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = {"released fixed model": "tab:blue", "our adapted/fine-tuned pretrained-CDS backbone": "tab:orange", "our adapted/fine-tuned generic pretrained backbone": "tab:green"}
    for kind, grp in df.groupby("kind"):
        ax.scatter(
            grp["intergenic_specificity"],
            grp["f1"],
            s=180 * np.clip(grp["predicted_gene_count_ratio_vs_reference"], 0.2, 2.5),
            c=cmap.get(kind, "tab:gray"),
            edgecolor="black",
            alpha=0.8,
            label=kind,
        )
        for _, row in grp.iterrows():
            ax.annotate(
                row["model_short"],
                (row["intergenic_specificity"], row["f1"]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=8,
            )

    ax.set_title("Clean Plant Aggregate: F1 vs Intergenic Specificity")
    ax.set_xlabel("Intergenic specificity")
    ax.set_ylabel("Constrained F1 at 0.01")
    ax.set_xlim(0.90, 1.00)
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3, linestyle="--")
    ax.axhline(0.75, color="grey", linestyle=":", linewidth=1)
    ax.axvline(0.99, color="red", linestyle="--", linewidth=1, label="specificity=0.99")
    ax.legend(loc="best", fontsize=8)

    fig.tight_layout()
    out = FIG_DIR / "result_clean_plant_frontier_scatter.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def plot_nt_unified_summary() -> None:
    df = parse_results_log_nt_points(RESULTS_MD)
    if df.empty:
        print("[warn] no NT-related rows parsed from results log")
        return

    df = df.copy().reset_index(drop=True)
    df["step"] = np.arange(len(df))
    df["f1_plot"] = pd.to_numeric(df["f1"], errors="coerce").clip(lower=0.0, upper=1.0)
    df["specificity_plot"] = pd.to_numeric(df["specificity"], errors="coerce").clip(lower=0.0, upper=1.0)
    df["fpr_plot"] = pd.to_numeric(df["fpr"], errors="coerce").clip(lower=0.0, upper=1.0)

    plot_df = df.dropna(subset=["specificity_plot", "f1_plot"]).copy().reset_index(drop=True)
    if plot_df.empty:
        print("[warn] parsed NT rows are missing both specificity and F1")
        return

    palette = {
        "NT-SEGNT": "#4E79A7",
        "NT-FRAGFIX": "#F28E2B",
        "M7-HELDOUT": "#59A14F",
        "M7-CLEAN": "#8CD17D",
        "M8-MULTICLASS": "#B07AA1",
        "M9-UNFREEZE": "#E15759",
        "M10-M9L12": "#B07C47",
        "GENERANNO": "#9C755F",
        "NT-OTHER": "#B0B0B0",
    }
    markers = {
        "NT-SEGNT": "o",
        "NT-FRAGFIX": "s",
        "M7-HELDOUT": "D",
        "M7-CLEAN": "P",
        "M8-MULTICLASS": "^",
        "M9-UNFREEZE": "v",
        "M10-M9L12": "X",
        "GENERANNO": "*",
        "NT-OTHER": "h",
    }

    # Build family summaries for bottom panels.
    family_summaries = []
    for family, sub in plot_df.groupby("family"):
        sub = sub.sort_values("step")
        if sub.empty:
            continue
        best_idx = sub["f1_plot"].idxmax()
        best = sub.loc[best_idx]
        latest = sub.iloc[-1]
        family_summaries.append(
            {
                "family": family,
                "best_f1": float(best["f1_plot"]),
                "best_specificity": float(best["specificity_plot"]),
                "best_fpr": float(best["fpr_plot"]) if not pd.isna(best["fpr_plot"]) else float("nan"),
                "latest_f1": float(latest["f1_plot"]),
                "latest_specificity": float(latest["specificity_plot"]),
                "latest_fpr": float(latest["fpr_plot"]) if not pd.isna(latest["fpr_plot"]) else float("nan"),
                "n": int(len(sub)),
            }
        )

    summary = pd.DataFrame(family_summaries)
    if summary.empty:
        print("[warn] NT rows parse succeeded but no family with usable points")
        return
    summary = summary.sort_values("best_f1", ascending=False).reset_index(drop=True)
    ordered_families = summary["family"].tolist()

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.2, 1.2], height_ratios=[2.2, 1.0])
    ax_tradeoff = fig.add_subplot(gs[0, :])
    ax_summary = fig.add_subplot(gs[1, 0])
    ax_fpr = fig.add_subplot(gs[1, 1])

    # Top panel: tradeoff trajectory (specificity vs F1), family-specific paths.
    for family in ordered_families:
        sub = plot_df[plot_df["family"] == family].sort_values("step")
        if sub.empty:
            continue
        color = palette.get(family, "tab:gray")
        marker = markers.get(family, "o")

        # point size encodes FPR (bigger means lower FPR); draw line for trajectory
        size = 60 + 180 * (1.0 - np.clip(sub["fpr_plot"].fillna(0.03), 0.0, 0.3))
        size = np.where(np.isfinite(size), size, 80)
        ax_tradeoff.plot(
            sub["specificity_plot"],
            sub["f1_plot"],
            marker=marker,
            linewidth=1.5,
            color=color,
            alpha=0.55,
            label=family,
            zorder=2,
        )
        ax_tradeoff.scatter(
            sub["specificity_plot"],
            sub["f1_plot"],
            s=size,
            facecolor=color,
            edgecolor="black",
            linewidth=0.3,
            alpha=0.85,
            marker=marker,
            zorder=3,
        )
        latest = sub.iloc[-1]
        ax_tradeoff.annotate(
            family,
            (latest["specificity_plot"], latest["f1_plot"]),
            xytext=(6, -12),
            textcoords="offset points",
            fontsize=8,
            color=color,
            alpha=0.8,
        )

    # top-3 F1 checkpoints as readable anchors
    for _, row in plot_df.nlargest(3, "f1_plot").iterrows():
        ax_tradeoff.annotate(
            row["run"][:16],
            (row["specificity_plot"], row["f1_plot"]),
            xytext=(4, 5),
            textcoords="offset points",
            fontsize=8,
            color="black",
        )

    ax_tradeoff.set_xlabel("intergenic specificity")
    ax_tradeoff.set_ylabel("Constrained/gene-body F1")
    ax_tradeoff.set_xlim(
        max(0.82, float(np.nanmin(plot_df["specificity_plot"])) - 0.008),
        min(1.0, float(np.nanmax(plot_df["specificity_plot"])) + 0.002),
    )
    ax_tradeoff.set_ylim(
        max(0.0, float(np.nanmin(plot_df["f1_plot"])) - 0.03),
        min(1.0, float(np.nanmax(plot_df["f1_plot"])) + 0.03),
    )
    ax_tradeoff.grid(alpha=0.2, linestyle="--")
    ax_tradeoff.axvline(0.99, color="#D62728", linestyle=":", linewidth=1, alpha=0.75, label="specificity=0.99")
    ax_tradeoff.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    handles, labels = ax_tradeoff.get_legend_handles_labels()
    uniq_map = {}
    for h, l in zip(handles, labels):
        if l and l not in uniq_map:
            uniq_map[l] = h
    ax_tradeoff.legend(uniq_map.values(), uniq_map.keys(), loc="best", fontsize=8)
    ax_tradeoff.set_title("NT trade-off map (specificity vs constrained/gene-body F1)")

    # Bottom-left: best checkpoint summary bars per family
    x = np.arange(len(summary))
    width = 0.35
    ax_summary.bar(
        x - width / 2,
        summary["best_f1"],
        width=width,
        color="#4C78A8",
        alpha=0.85,
        label="best F1",
    )
    ax_summary.bar(
        x + width / 2,
        summary["best_specificity"],
        width=width,
        color="#F28E2B",
        alpha=0.85,
        label="specificity@bestF1",
    )
    ax_summary.set_title("Best checkpoint comparison by family")
    ax_summary.set_ylabel("metric")
    ax_summary.set_ylim(0.0, 1.05)
    ax_summary.set_xticks(x)
    ax_summary.set_xticklabels(summary["family"].tolist(), rotation=28, ha="right", fontsize=8)
    ax_summary.grid(alpha=0.25, axis="y", linestyle="--")
    ax_summary.legend(fontsize=8, loc="best")
    for idx, row in summary.iterrows():
        ax_summary.text(idx - width / 2, row["best_f1"] + 0.013, f"{row['best_f1']:.3f}", fontsize=7, ha="center")
        ax_summary.text(idx + width / 2, row["best_specificity"] + 0.013, f"{row['best_specificity']:.3f}", fontsize=7, ha="center")

    # Bottom-right: FPR spread by family with hard lines for the guardrail thresholds.
    fpr_groups = []
    fpr_labels = []
    fpr_max = 0.0
    for family in ordered_families:
        vals = plot_df[plot_df["family"] == family]["fpr_plot"].dropna().to_numpy()
        if len(vals) == 0:
            continue
        fpr_groups.append(vals)
        fpr_labels.append(family)
        fpr_max = max(fpr_max, float(np.nanmax(vals)))

    if fpr_groups:
        bp = ax_fpr.boxplot(
            fpr_groups,
            widths=0.6,
            patch_artist=True,
            showmeans=True,
            positions=np.arange(len(fpr_groups)),
        )
        ax_fpr.set_xticks(np.arange(len(fpr_groups)))
        ax_fpr.set_xticklabels(fpr_labels, rotation=28, ha="right", fontsize=8)
        for box, family in zip(bp["boxes"], fpr_labels):
            box.set_facecolor(palette.get(family, "#B0B0B0"))
            box.set_alpha(0.35)
        ax_fpr.set_ylim(0.0, max(0.06, fpr_max * 1.15 + 1e-6))
        ax_fpr.set_ylabel("intergenic FPR")
    else:
        ax_fpr.text(0.5, 0.5, "No finite FPR parsed for NT points", ha="center", va="center", transform=ax_fpr.transAxes)
        ax_fpr.set_ylabel("intergenic FPR")
        ax_fpr.set_ylim(0.0, 0.05)

    ax_fpr.set_title("FPR distribution by family (lower is better)")
    ax_fpr.grid(alpha=0.2, linestyle="--", axis="y")
    line_002 = ax_fpr.axhline(
        0.02,
        color="#D62728",
        linestyle="--",
        linewidth=1,
        label="FPR=0.02",
    )
    line_001 = ax_fpr.axhline(
        0.01,
        color="#9467BD",
        linestyle=":",
        linewidth=1,
        label="FPR=0.01",
    )
    if fpr_groups:
        ax_fpr.legend(handles=[line_002, line_001], fontsize=7, loc="best")

    fig.tight_layout()
    out = FIG_DIR / "result_nt_unified_summary.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def plot_clean_plant_per_species() -> None:
    per_species_path = (
        REPORTS_DIR / "M19-COMPARABILITY-EVIDENCE" / "clean_plant_per_species.csv"
    )
    if not per_species_path.exists():
        print(f"[warn] missing {per_species_path}, skip per-species plot")
        return

    df = pd.read_csv(per_species_path)
    df["model_short"] = df["model"].str.split(",").str[0]
    model_order = df["model_short"].unique().tolist()
    species = df["species"].unique().tolist()

    bar_w = 0.18
    x = np.arange(len(model_order))
    fig, (ax_f1, ax_fpr) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    for si, sp in enumerate(species):
        sub = df[df["species"] == sp]
        offsets = x + si * bar_w - bar_w * (len(species) - 1) / 2
        ax_f1.bar(
            offsets,
            sub.set_index("model_short").reindex(model_order)["gene_body_F1_unconstrained"],
            width=bar_w,
            label=sp,
            alpha=0.9,
        )
        ax_fpr.bar(
            offsets,
            sub.set_index("model_short").reindex(model_order)["intergenic_FPR"],
            width=bar_w,
            label=sp,
            alpha=0.9,
        )

    ax_f1.set_ylabel("Gene-body F1")
    ax_f1.set_ylim(0, 1.0)
    ax_f1.set_title("Clean Plant Per-Species F1")
    ax_f1.grid(alpha=0.2, axis="y")

    ax_fpr.set_ylabel("Intergenic FPR")
    ax_fpr.set_title("Clean Plant Per-Species Intergenic FPR")
    ax_fpr.set_xlabel("Model")
    ax_fpr.set_xticks(x)
    ax_fpr.set_xticklabels(model_order, rotation=22, ha="right")
    ax_fpr.grid(alpha=0.2, axis="y")
    ax_fpr.set_ylim(0, 0.03)

    ax_f1.legend(title="Species", bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=8)
    ax_fpr.legend(title="Species", bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=8)

    fig.tight_layout()
    out = FIG_DIR / "result_clean_plant_per_species.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def plot_results_log_timeline() -> None:
    df = parse_results_log_aggregates(RESULTS_MD)
    if df.empty:
        print("[warn] no parseable aggregate metrics found in results log")
        return

    df = df.reset_index(drop=True)
    fig, (ax_f1, ax_spec) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    ax_f1.plot(df["experiment"], df["f1"], marker="o")
    ax_f1.set_title("Results log aggregate constrained F1 trajectory")
    ax_f1.set_ylabel("Constrained F1 (best available)")
    ax_f1.set_ylim(0, 1.0)
    ax_f1.grid(alpha=0.3, linestyle="--")

    ax_spec.plot(df["experiment"], df["specificity"], marker="o", color="tab:orange")
    ax_spec.set_title("Results log aggregate intergenic specificity trajectory")
    ax_spec.set_ylabel("Intergenic specificity")
    ax_spec.set_xlabel("Experiment")
    ax_spec.set_ylim(0.9, 1.0)
    ax_spec.grid(alpha=0.3, linestyle="--")
    ax_spec.set_xticks(range(len(df["experiment"])))
    ax_spec.set_xticklabels(df["experiment"], rotation=55, ha="right")

    fig.tight_layout()
    out = FIG_DIR / "result_log_trajectory.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def plot_oracle_tradeoff() -> None:
    metric_files = sorted(
        (REPORTS_DIR / "M18-MULTICLADE-ORACLE-CALIBRATION" / "gallus_gallus" / "metrics").glob(
            "*.metrics.json"
        )
    )
    if not metric_files:
        print("[warn] missing oracle metric files, skip oracle tradeoff plot")
        return

    rows = []
    for path in metric_files:
        payload = json.loads(path.read_text(encoding="utf-8"))
        rows.append(
            {
                "decode": path.stem.replace(".metrics", ""),
                "gene_body_F1": payload.get("gene_body_F1_unconstrained", float("nan")),
                "constrained_F1": payload.get("constrained_gene_body_F1", float("nan")),
                "specificity": payload.get("intergenic_specificity", float("nan")),
                "fpr": payload.get("intergenic_FPR", float("nan")),
                "gene_count_ratio": payload.get("predicted_gene_count_ratio_vs_reference", float("nan")),
                "predicted_genes": payload.get("predicted_gene_count", float("nan")),
                "reference_genes": payload.get("reference_gene_count", float("nan")),
            }
        )
    oracle = pd.DataFrame(rows)
    oracle = oracle.replace([np.inf, -np.inf], np.nan).dropna(subset=["gene_body_F1", "fpr"])

    # Keep all points; use small marker for dense cloud.
    fig, ax = plt.subplots(figsize=(10, 6))
    sc = ax.scatter(
        oracle["fpr"],
        oracle["gene_body_F1"],
        c=oracle["gene_count_ratio"],
        cmap="viridis",
        s=35,
        alpha=0.8,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Predicted gene count ratio vs reference")
    ax.axvline(0.02, color="red", linestyle="--", linewidth=1, label="FPR=0.02")
    ax.axvline(0.01, color="orange", linestyle=":", linewidth=1, label="FPR=0.01")
    ax.set_title("M18 Oracle Calibration Sweep (Gallus)")
    ax.set_xlabel("Intergenic FPR")
    ax.set_ylabel("Unconstrained gene-body F1")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(loc="best")

    # Mark one good and one bad frontier point by constrained F1 threshold.
    if not oracle.empty:
        idx_best = oracle["gene_body_F1"].idxmax()
        idx_worst = oracle["fpr"].idxmax()
        ax.scatter(
            [oracle.loc[idx_best, "fpr"]],
            [oracle.loc[idx_best, "gene_body_F1"]],
            marker="*",
            s=180,
            edgecolor="black",
            c="gold",
            label="Best F1",
        )
        ax.scatter(
            [oracle.loc[idx_worst, "fpr"]],
            [oracle.loc[idx_worst, "gene_body_F1"]],
            marker="X",
            s=160,
            edgecolor="black",
            c="red",
            label="Worst FPR",
        )

    fig.tight_layout()
    out = FIG_DIR / "result_oracle_gallus_tradeoff.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def plot_m13_seed_diagnostics() -> None:
    tsv = REPORTS_DIR / "M13_FAILURE_SANITY" / "per_seed.tsv"
    if not tsv.exists():
        print(f"[warn] missing {tsv}, skip M13 seed diagnostics")
        return

    df = pd.read_csv(tsv, sep="\t")
    df["run"] = df["family"] + "-s" + df["seed"].astype(str)
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    ax[0].plot(df["run"], df["val_gbf1"], marker="o", label="val gbF1")
    ax[0].plot(df["run"], df["test_gbf1"], marker="s", label="test gbF1")
    ax[0].set_ylabel("gene-body gbF1")
    ax[0].set_ylim(0, 1)
    ax[0].grid(alpha=0.3, linestyle="--")
    ax[0].legend(loc="best")

    ax[1].plot(df["run"], df["val_fpr"], marker="o", label="val FPR")
    ax[1].plot(df["run"], df["test_fpr"], marker="s", label="test FPR")
    ax[1].set_ylabel("FPR")
    ax[1].set_xlabel("run")
    ax[1].set_ylim(0, 0.1)
    ax[1].grid(alpha=0.3, linestyle="--")
    ax[1].legend(loc="best")
    ax[1].set_xticks(range(len(df)))
    ax[1].set_xticklabels(df["run"].tolist(), rotation=30, ha="right")

    fig.suptitle("M13 Failure-Sanity seeds (raw diagnostics)")
    fig.tight_layout()
    out = FIG_DIR / "result_m13_seed_diagnostics.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[ok] wrote {out}")


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    plot_clean_plant_frontier()
    plot_clean_plant_per_species()
    plot_results_log_timeline()
    plot_oracle_tradeoff()
    plot_nt_unified_summary()
    plot_m13_seed_diagnostics()

    print("[done] all requested figures generated")


if __name__ == "__main__":
    main()
