"""TB-GBF1-MULTICLASS-M8 CK4 collection: mc-candidate(5) vs 3c-candidate(5) vs anchor(3) on the
CLEAN held-out plants {arabidopsis, rice}. KEY question: does multi-class RECOVER gbF1 (> 3c and
toward/over the raw-DNA anchor) WHILE holding intergenic_specificity >=~0.95? base-weighted mean±std
+ per-species. Run: python scripts/_collect_m8.py"""
import json, statistics as st

KEYS = ("intergenic_specificity", "macro_intergenic_specificity",
        "gene_body_F1_unconstrained", "predicted_gene_count_ratio_vs_reference")
SHORT = {"intergenic_specificity": "spec", "macro_intergenic_specificity": "macro",
         "gene_body_F1_unconstrained": "gbF1", "predicted_gene_count_ratio_vs_reference": "gcount"}


def load(e):
    try:
        return json.load(open(f"outputs/{e}/metrics/metrics.json"))
    except Exception:
        return None


def agg(exps):
    rows = [r for r in (load(e) for e in exps) if r]
    out = {}
    for k in KEYS:
        v = [r[k] for r in rows if r.get(k) is not None]
        out[k] = (round(st.mean(v), 4), round(st.pstdev(v), 4) if len(v) > 1 else 0.0, len(v))
    return out, rows


def persp(rows, metric):
    acc = {}
    for r in rows:
        for sp, v in r.get("per_species", {}).items():
            key = "arab" if "arabidopsis" in sp else ("rice" if "oryza" in sp else sp)
            acc.setdefault(key, []).append(v.get(metric))
    return {k: round(st.mean([x for x in vs if x is not None]), 4) for k, vs in acc.items()}


groups = {
    "anchor(raw-DNA,3c)": [f"SCREENREF-tiberius_like-m8clean-s{s}" for s in range(3)],
    "3c-candidate(M7 cfg)": [f"M8-3C-CAND-s{s}" for s in range(5)],
    "mc-candidate(M8 multi)": [f"M8-MC-CAND-s{s}" for s in range(5)],
}
res = {}
for name, exps in groups.items():
    a, rows = agg(exps)
    res[name] = a
    print(f"\n=== {name} (n={a['intergenic_specificity'][2]}) ===")
    for k in KEYS:
        print(f"  {SHORT[k]:7s} {a[k][0]} +-{a[k][1]}")
    print("  per-sp spec:", persp(rows, "intergenic_specificity"), "| per-sp gbF1:", persp(rows, "gene_body_F1_unconstrained"))

print("\n=== KEY M8 VERDICT (clean held-out {arabidopsis,rice}) ===")
A = res["anchor(raw-DNA,3c)"]; C3 = res["3c-candidate(M7 cfg)"]; MC = res["mc-candidate(M8 multi)"]
ag, c3g, mcg = A["gene_body_F1_unconstrained"][0], C3["gene_body_F1_unconstrained"][0], MC["gene_body_F1_unconstrained"][0]
asp, c3sp, mcsp = A["intergenic_specificity"][0], C3["intergenic_specificity"][0], MC["intergenic_specificity"][0]
print(f"  gbF1:  anchor {ag} | 3c {c3g} | mc {mcg}")
print(f"  spec:  anchor {asp} | 3c {c3sp} | mc {mcsp}")
print(f"  -> mc gbF1 > 3c gbF1 (recovery)? {'YES' if mcg > c3g else 'NO'} ({mcg-c3g:+.4f})")
print(f"  -> mc gbF1 >= anchor gbF1 (bar)? {'YES' if mcg >= ag else 'NO'} ({mcg-ag:+.4f})")
print(f"  -> mc spec held >=0.95? {'YES' if mcsp >= 0.95 else 'NO'} ({mcsp})")
print(f"  -> mc gcount in [1.0,1.25]? {MC['predicted_gene_count_ratio_vs_reference'][0]}")
