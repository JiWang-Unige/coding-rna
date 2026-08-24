import json, numpy as np, os
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
seeds=[0,1,2]
rows=[]
for s in seeds:
    p=f"{ROOT}/outputs/FP-SEGMENTNT-PROBE-M1-convlstm-s{s}/metrics/metrics.json"
    rows.append(json.load(open(p)))
def mean(k): return float(np.mean([r[k] for r in rows]))
def std(k): return float(np.std([r[k] for r in rows]))
print("=== seed-mean (base-weighted), 3 seeds ===")
for k in ["intergenic_specificity","macro_intergenic_specificity","gene_body_F1_unconstrained",
          "gene_body_precision","gene_body_recall","intergenic_FPR","predicted_gene_count_ratio_vs_reference"]:
    vals=[round(r[k],4) for r in rows]
    print("  %-42s mean=%.4f std=%.4f per-seed=%s" % (k, mean(k), std(k), vals))
print("=== per-species per-seed ===")
for s,r in zip(seeds,rows):
    ps=r.get("per_species",{})
    print("  seed%d:" % s)
    for sp,v in ps.items():
        print("    %-30s spec=%.4f gbF1=%.4f FPR=%.4f gcount=%s/%s" % (
            sp, v.get("intergenic_specificity") or 0, v.get("gene_body_F1_unconstrained") or 0,
            v.get("intergenic_FPR") or 0, v.get("predicted_gene_count"), v.get("reference_gene_count")))
