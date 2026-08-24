#!/usr/bin/env python3
"""One-off: recompute the same-budget screen ladder under the NEW intergenic ruler
(revise-goal 2026-06-11: intergenic = complement of FULL-transcript span incl UTR;
PRIMARY metric = intergenic_specificity = 1 - intergenic_FPR).

Re-evaluates EXISTING prediction GFFs (no GPU, no retraining) against each run's
own held-out eval_subsets, aggregates base-weighted + macro across the 2 pilot
species, and prints the ladder so screen_anchor can be recalibrated.
"""
import json
import subprocess
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
SPECIES = ["saccharomyces_cerevisiae", "drosophila_melanogaster"]

# (label, run-dir glob over seeds)
RUNS = {
    "FLOOR(ORF)":            ["FLOOR-SCREEN-M1"],
    "tiberius_like(anchor)": [f"SCREENREF-tiberius_like-s{s}" for s in (0, 1, 2)],
    "helixer_like":          [f"SCREENREF-helixer_like-s{s}" for s in (0, 1, 2)],
    "CONSTR(post-proc)":     [f"SCREENREF-tiberius_like-constrained-s{s}" for s in (0, 1, 2)],
    "CRF-vec(learned)":      [f"SCREENREF-tiberius_like-crf-s{s}" for s in (0, 1, 2)],
}


def eval_run(run_dir: Path, workdir: Path):
    """Re-eval both species for one run, aggregate, return the aggregate dict (or None)."""
    sp_jsons = []
    for sp in SPECIES:
        ref = run_dir / "eval_subsets" / sp / "reference.gff3"
        pred = run_dir / "predictions" / sp / "..." if False else run_dir / "predictions" / f"{sp}.gff"
        genome = run_dir / "eval_subsets" / sp / "genome.fa"
        if not (ref.exists() and pred.exists() and genome.exists()):
            return None
        out = workdir / f"{run_dir.name}_{sp}.json"
        subprocess.run([PY, str(ROOT / "scripts/eval_gene_body_mask.py"),
                        "--reference-gtf", str(ref), "--prediction-gtf", str(pred),
                        "--genome-fasta", str(genome), "--output-json", str(out),
                        "--experiment-id", f"REANCHOR_{run_dir.name}_{sp}",
                        "--profile", "screen", "--span-mode", "cds"], check=True)
        sp_jsons.append(str(out))
    agg = workdir / f"{run_dir.name}_AGG.json"
    subprocess.run([PY, str(ROOT / "scripts/aggregate_gene_body_metrics.py"),
                    "--metrics", *sp_jsons, "--output-json", str(agg),
                    "--experiment-id", f"REANCHOR_{run_dir.name}", "--profile", "screen"], check=True)
    return json.load(open(agg))


def main():
    workdir = Path("/tmp/reanchor_newruler")
    workdir.mkdir(parents=True, exist_ok=True)
    print(f"{'run':<24}{'seeds':>6}{'spec_bw':>9}{'spec_macro':>11}{'FPR_bw':>8}"
          f"{'gbF1_bw':>9}{'gcount_ratio':>13}")
    print("-" * 80)
    ladder = {}
    for label, run_names in RUNS.items():
        specs_bw, specs_macro, fprs, f1s, ratios, n = [], [], [], [], [], 0
        for rn in run_names:
            d = eval_run(ROOT / "outputs" / rn, workdir)
            if d is None:
                continue
            n += 1
            specs_bw.append(d["intergenic_specificity"])
            if d.get("macro_intergenic_specificity") is not None:
                specs_macro.append(d["macro_intergenic_specificity"])
            fprs.append(d["intergenic_FPR"])
            f1s.append(d["gene_body_F1_unconstrained"])
            ratios.append(d["predicted_gene_count_ratio_vs_reference"])
        if not n:
            print(f"{label:<24}{'MISSING':>6}")
            continue
        row = dict(seeds=n, spec_bw=mean(specs_bw),
                   spec_macro=mean(specs_macro) if specs_macro else float("nan"),
                   fpr_bw=mean(fprs), gbf1_bw=mean(f1s), gcount=mean(ratios),
                   spec_bw_perseed=specs_bw)
        ladder[label] = row
        print(f"{label:<24}{n:>6}{row['spec_bw']:>9.4f}{row['spec_macro']:>11.4f}"
              f"{row['fpr_bw']:>8.4f}{row['gbf1_bw']:>9.4f}{row['gcount']:>13.2f}")
    print("-" * 80)
    json.dump(ladder, open(workdir / "ladder_summary.json", "w"), indent=2)
    print(f"saved: {workdir / 'ladder_summary.json'}")


if __name__ == "__main__":
    main()
