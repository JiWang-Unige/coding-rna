#!/usr/bin/env python3
"""M9-DEEP: does deeper unfreeze push FPR<=0.02 (spec>=0.98) to unlock constrained_gbF1?
Arms L6/L8/L12 vs prior L4 baseline (gbF1 0.8759 spec 0.9754 FPR 0.0246).
Barrier: screen constrained threshold FPR<=0.02. ANNEVO ceiling gbF1 0.898.
"""
import json, os, sys
ROOT = "/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ARMS = {6: "unfreeze-L6", 8: "unfreeze-L8", 12: "unfreeze-L12"}
L4 = dict(spec=0.9754, gbF1=0.8759, fpr=0.0246, cgbF1=0.0, gcount=0.820)
def load(n):
    p = f"{ROOT}/outputs/M9-UNFREEZE-L{n}-s0/metrics/metrics.json"
    return json.load(open(p)) if os.path.exists(p) else None
rows = {}
for n in ARMS:
    m = load(n)
    if m is None:
        print(f"[WARN] L{n}-s0 metrics absent (job not done)"); continue
    rows[n] = dict(spec=m.get("intergenic_specificity", -1), gbF1=m.get("gene_body_F1_unconstrained", -1),
                   cgbF1=m.get("constrained_gene_body_F1", -1), fpr=m.get("intergenic_FPR", -1),
                   gcount=m.get("predicted_gene_count_ratio_vs_reference", -1))
if not rows:
    print("No deep-arm metrics yet."); sys.exit(0)
print("\n=== M9-DEEP unfreeze depth sweep (arabidopsis, base-w, 1 seed) ===")
print(f"{'arm':<16}{'spec':>8}{'gbF1':>8}{'cgbF1':>8}{'FPR':>8}{'gcount':>8}{'FPR<=.02?':>10}")
print(f"{'L4 (prior)':<16}{L4['spec']:>8.4f}{L4['gbF1']:>8.4f}{L4['cgbF1']:>8.4f}{L4['fpr']:>8.4f}{L4['gcount']:>8.3f}{'NO':>10}")
broke = []
for n in ARMS:
    if n not in rows: print(f"{ARMS[n]:<16}{'(pending)':>8}"); continue
    r = rows[n]; ok = r['fpr'] <= 0.02
    if ok: broke.append(n)
    print(f"{ARMS[n]:<16}{r['spec']:>8.4f}{r['gbF1']:>8.4f}{r['cgbF1']:>8.4f}{r['fpr']:>8.4f}{r['gcount']:>8.3f}{('YES' if ok else 'no'):>10}")
print(f"\nANNEVO ceiling gbF1=0.898 | barrier: FPR<=0.02 unlocks constrained_gbF1")
print("\n=== VERDICT ===")
if broke:
    best = max(broke, key=lambda n: rows[n]['cgbF1'])
    print(f"BARRIER BROKEN by L{broke}: deeper unfreeze pushed FPR<=0.02. Best constrained_gbF1=L{best} {rows[best]['cgbF1']:.4f}. -> validate_goal should now PASS; promote to multi-seed+cross-species CI.")
else:
    print("BARRIER HELD: even L6/L8/L12 FPR still >0.02. Need full-unfreeze OR FP-aware objective tweak (higher fp-lambda) OR constrained-decode threshold pass. Report monotonic FPR trend to decide.")
