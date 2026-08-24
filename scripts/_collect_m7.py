"""REANCHOR-HELDOUT-M7 CK6 collection + Pareto judgment.
Collect candidate FP-FRAGFIX-CONSTR-ho s0-4 (5 seeds) + anchor SCREENREF-tiberius_like-ho s0-2
(3 seeds) + ANNEVO ceiling. Report base-weighted mean+-std for spec/macro/gbF1/gene_count + per-
species spec. Pareto verdict vs held-out anchor. Run: python scripts/_collect_m7.py
"""
import json, statistics as st

KEYS = ("intergenic_specificity", "macro_intergenic_specificity",
        "gene_body_F1_unconstrained", "predicted_gene_count_ratio_vs_reference")


def load(exp):
    try:
        return json.load(open(f"outputs/{exp}/metrics/metrics.json"))
    except Exception:
        return None


def agg(exps):
    rows = [load(e) for e in exps]
    rows = [r for r in rows if r]
    out = {}
    for k in KEYS:
        vals = [r.get(k) for r in rows if r.get(k) is not None]
        out[k] = (round(st.mean(vals), 4), round(st.pstdev(vals), 4) if len(vals) > 1 else 0.0,
                  [round(v, 4) for v in vals])
    return out, rows


def persp(rows, metric="intergenic_specificity"):
    acc = {}
    for r in rows:
        for sp, v in r.get("per_species", {}).items():
            short = sp.split("_")[0] if "_" in sp else sp
            for k in ("arabidopsis", "gallus"):
                if k in sp:
                    short = k
            acc.setdefault(short, []).append(v.get(metric))
    return {k: round(st.mean([x for x in vs if x is not None]), 4) for k, vs in acc.items()}


cand_exps = [f"FP-FRAGFIX-CONSTR-ho-s{s}" for s in range(5)]
anch_exps = [f"SCREENREF-tiberius_like-ho-s{s}" for s in range(3)]
cand, cand_rows = agg(cand_exps)
anch, anch_rows = agg(anch_exps)
ceil = load("REANCHOR-CEILING-ANNEVO-M7")

print("=== held-out anchor (tiberius_like, 3 seeds, base-w) ===")
for k in KEYS:
    print(f"  {k}: {anch[k][0]} +-{anch[k][1]}  seeds={anch[k][2]}")
print("  per-species spec:", persp(anch_rows))
print("\n=== candidate FP-FRAGFIX-CONSTR (5 seeds, base-w) ===")
for k in KEYS:
    print(f"  {k}: {cand[k][0]} +-{cand[k][1]}  seeds={cand[k][2]}")
print("  per-species spec:", persp(cand_rows))
print("  per-species gbF1:", persp(cand_rows, "gene_body_F1_unconstrained"))
print("  per-species gcount:", persp(cand_rows, "predicted_gene_count_ratio_vs_reference"))
if ceil:
    print("\n=== ANNEVO ceiling (held-out, test-chroms) ===")
    print("  spec", round(ceil.get("intergenic_specificity", -1), 4),
          "macro", round(ceil.get("macro_intergenic_specificity", -1), 4),
          "gbF1", round(ceil.get("gene_body_F1_unconstrained", -1), 4))

# Pareto verdict
A_spec, C_spec = anch["intergenic_specificity"][0], cand["intergenic_specificity"][0]
A_macro = anch["macro_intergenic_specificity"][0]
C_macro = cand["macro_intergenic_specificity"][0]
C_gbf1 = cand["gene_body_F1_unconstrained"][0]
C_gc = cand["predicted_gene_count_ratio_vs_reference"][0]
FLOOR = 0.5276
print("\n=== PARETO VERDICT (held-out) ===")
print(f"  spec  : cand {C_spec} {'>' if C_spec>A_spec else '<='} anchor {A_spec}  -> {'PASS' if C_spec>A_spec else 'FAIL'}")
print(f"  macro : cand {C_macro} {'>=' if C_macro>=A_macro else '<'} anchor-macro-gate {A_macro}  -> {'PASS' if C_macro>=A_macro else 'FAIL'}")
print(f"  gbF1  : cand {C_gbf1} {'>=' if C_gbf1>=FLOOR else '<'} floor {FLOOR}  -> {'PASS' if C_gbf1>=FLOOR else 'FAIL'}")
print(f"  gcount: cand {C_gc} in [1.0,1.25]? -> {'PASS' if 1.0<=C_gc<=1.25 else 'WRINKLE (band re-select on VAL)'}")
