"""Collect TA-FOUNDATION-DECODER-M4 batch: per-candidate seed-mean +- std + per-species, vs anchor.
Run: python3 scripts/_collect_m4.py"""
import json, os, numpy as np
ROOT = "/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
CANDS = ["FPLOSS", "FUSION", "CRF"]
SEEDS = [0, 1, 2, 3, 4]
ANCHOR_SPEC, ANCHOR_MACRO, F1_FLOOR, MACRO_GATE = 0.8710, 0.8278, 0.5276, 0.7978
KEYS = ["intergenic_specificity", "macro_intergenic_specificity", "gene_body_F1_unconstrained",
        "intergenic_FPR", "predicted_gene_count_ratio_vs_reference"]


def load(cand, s):
    p = f"{ROOT}/outputs/FP-SEGNT-{cand}-s{s}/metrics/metrics.json"
    return json.load(open(p)) if os.path.exists(p) else None


print(f"{'cand':<8}{'n':>3}{'spec_bw(±)':>18}{'macro':>9}{'gbF1':>8}{'FPR':>8}{'gcount':>8}{'PASS?':>8}")
print("-" * 72)
for cand in CANDS:
    rows = [load(cand, s) for s in SEEDS]
    rows = [r for r in rows if r]
    if not rows:
        print(f"{cand:<8}  MISSING"); continue
    m = {k: float(np.mean([r[k] for r in rows])) for k in KEYS}
    sd = {k: float(np.std([r[k] for r in rows])) for k in KEYS}
    passes = (m["intergenic_specificity"] > ANCHOR_SPEC and
              m["gene_body_F1_unconstrained"] >= F1_FLOOR and
              m["macro_intergenic_specificity"] >= MACRO_GATE)
    print(f"{cand:<8}{len(rows):>3}{m['intergenic_specificity']:>10.4f}(±{sd['intergenic_specificity']:.3f})"
          f"{m['macro_intergenic_specificity']:>9.4f}{m['gene_body_F1_unconstrained']:>8.4f}"
          f"{m['intergenic_FPR']:>8.4f}{m['predicted_gene_count_ratio_vs_reference']:>8.2f}"
          f"{('YES' if passes else 'no'):>8}")
    print(f"         per-seed spec: {[round(r['intergenic_specificity'],4) for r in rows]}")
    # per-species (from per_species of each seed, averaged)
    for sp in ("saccharomyces_cerevisiae", "drosophila_melanogaster"):
        sv = [r.get("per_species", {}).get(sp, {}) for r in rows]
        sv = [v for v in sv if v]
        if sv:
            spec = np.mean([v.get("intergenic_specificity") or 0 for v in sv])
            f1 = np.mean([v.get("gene_body_F1_unconstrained") or 0 for v in sv])
            print(f"           {sp[:12]:<12} spec={spec:.4f} gbF1={f1:.4f}")
print("-" * 72)
print(f"anchor: spec_bw {ANCHOR_SPEC} / macro {ANCHOR_MACRO} / F1_floor {F1_FLOOR} / macro_gate {MACRO_GATE} / ceiling 0.9917")
print("PASS = spec_bw>anchor AND gbF1>=floor AND macro>=gate (screen direction-selection; NON-CLAIM)")
