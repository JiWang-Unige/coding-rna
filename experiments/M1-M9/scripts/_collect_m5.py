"""Collect TA-COHERENCE-FIX-M5: re-eval the 5-seed anchor (tiberius_like s0-4) under the NEW ruler,
collect FP-FRAGFIX-CONSTR 5 seeds, paired comparison + fragmentation. python3 scripts/_collect_m5.py"""
import json, os, subprocess, sys, numpy as np
ROOT = "/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
PY = sys.executable
SP = ["saccharomyces_cerevisiae", "drosophila_melanogaster"]
W = "/tmp/m5_collect"; os.makedirs(W, exist_ok=True)
F1_FLOOR, MACRO_GATE = 0.5276, 0.7978


def neweval(run):
    """Re-eval a run's predictions under the NEW ruler -> aggregate dict (or None)."""
    d = f"{ROOT}/outputs/{run}"
    sj = []
    for sp in SP:
        ref = f"{d}/eval_subsets/{sp}/reference.gff3"; gen = f"{d}/eval_subsets/{sp}/genome.fa"
        pred = f"{d}/predictions/{sp}.gff"
        if not all(os.path.exists(x) for x in (ref, gen, pred)):
            return None
        o = f"{W}/{run}_{sp}.json"
        subprocess.run([PY, f"{ROOT}/scripts/eval_gene_body_mask.py", "--reference-gtf", ref,
                        "--prediction-gtf", pred, "--genome-fasta", gen, "--output-json", o,
                        "--experiment-id", f"{run}_{sp}", "--profile", "screen", "--span-mode", "cds"], check=True)
        sj.append(o)
    agg = f"{W}/{run}_AGG.json"
    subprocess.run([PY, f"{ROOT}/scripts/aggregate_gene_body_metrics.py", "--metrics", *sj,
                    "--output-json", agg, "--experiment-id", run, "--profile", "screen"], check=True)
    return json.load(open(agg))


def collect(runs):
    specs, macros, f1s, gcs = [], [], [], []
    for r in runs:
        m = neweval(r) if r.startswith("SCREENREF") else (
            json.load(open(f"{ROOT}/outputs/{r}/metrics/metrics.json"))
            if os.path.exists(f"{ROOT}/outputs/{r}/metrics/metrics.json") else None)
        if not m:
            continue
        specs.append(m["intergenic_specificity"]); macros.append(m.get("macro_intergenic_specificity") or 0)
        f1s.append(m["gene_body_F1_unconstrained"]); gcs.append(m["predicted_gene_count_ratio_vs_reference"])
    return specs, macros, f1s, gcs


anchor_runs = [f"SCREENREF-tiberius_like-s{s}" for s in range(5)]
constr_runs = [f"FP-FRAGFIX-CONSTR-s{s}" for s in range(5)]
a_spec, a_macro, a_f1, a_gc = collect(anchor_runs)
c_spec, c_macro, c_f1, c_gc = collect(constr_runs)

print("=== 5-seed ANCHOR (tiberius_like, NEW ruler re-eval) ===")
print(f"  spec per-seed {[round(x,4) for x in a_spec]} mean={np.mean(a_spec):.4f}±{np.std(a_spec):.3f} | "
      f"macro {np.mean(a_macro):.4f} | gbF1 {np.mean(a_f1):.4f} | gcount {np.mean(a_gc):.2f}")
print("=== FP-FRAGFIX-CONSTR (FPLOSS + constrained post-proc) ===")
print(f"  spec per-seed {[round(x,4) for x in c_spec]} mean={np.mean(c_spec):.4f}±{np.std(c_spec):.3f} | "
      f"macro {np.mean(c_macro):.4f} | gbF1 {np.mean(c_f1):.4f} | gcount {np.mean(c_gc):.2f}")
if a_spec and c_spec:
    anc_mean = np.mean(a_spec)
    paired = None
    if len(a_spec) == len(c_spec):
        diff = np.array(c_spec) - np.array(a_spec)
        paired = (float(diff.mean()), float(diff.std()))
    cm, cM, cf, cg = np.mean(c_spec), np.mean(c_macro), np.mean(c_f1), np.mean(c_gc)
    print("=== GATE (CONSTR vs 5-seed anchor) ===")
    print(f"  spec {cm:.4f} > anchor {anc_mean:.4f} ? {cm > anc_mean}")
    print(f"  gbF1 {cf:.4f} >= {F1_FLOOR} ? {cf >= F1_FLOOR}")
    print(f"  macro {cM:.4f} >= {MACRO_GATE} ? {cM >= MACRO_GATE}")
    print(f"  gene_count_ratio {cg:.2f} -> toward <=1.25 ? {cg <= 1.25} (FPLOSS-no-postproc was 2.25)")
    if paired:
        print(f"  paired diff (CONSTR - anchor) mean={paired[0]:.4f} ± {paired[1]:.3f}")
    print(f"  PROMOTE-READY ? {cm > anc_mean and cf >= F1_FLOOR and cM >= MACRO_GATE and cg <= 1.25}")
