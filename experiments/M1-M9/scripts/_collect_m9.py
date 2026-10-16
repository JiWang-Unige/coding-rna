#!/usr/bin/env python3
"""M9-CK4/CK5: collect 3-arm unfreeze results + judge primary_progress_gate.
3 arms (1 seed, single-species arabidopsis): L0=frozen CONTROL, L2/L4=unfreeze.
Gate: unfreeze gbF1 > frozen-L0 control AND > 14-elem 3c 0.7392 (M8) AND spec>=0.93.
Anchors: M8 3c gbF1 0.7392 / spec 0.966 ; ANNEVO ceiling gbF1 0.898.
"""
import json, os, sys
ROOT = "/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
ARMS = {0: "frozen-L0 (CONTROL)", 2: "unfreeze-L2", 4: "unfreeze-L4"}
M8_3C_GBF1, M8_3C_SPEC, ANNEVO_CEIL = 0.7392, 0.966, 0.898

def load(n):
    p = f"{ROOT}/outputs/M9-UNFREEZE-L{n}-s0/metrics/metrics.json"
    return json.load(open(p)) if os.path.exists(p) else None

rows = {}
for n in ARMS:
    m = load(n)
    if m is None:
        print(f"[WARN] L{n}-s0 metrics absent (job not done)")
        continue
    rows[n] = dict(
        spec=m.get("intergenic_specificity", -1),
        gbF1=m.get("gene_body_F1_unconstrained", -1),
        cgbF1=m.get("constrained_gene_body_F1", -1),
        fpr=m.get("intergenic_FPR", -1),
        gcount=m.get("predicted_gene_count_ratio_vs_reference", -1),
    )

if not rows:
    print("No arm metrics yet.")
    sys.exit(0)

print("\n=== M9 CK4 3-arm (arabidopsis, base-w, 1 seed) ===")
print(f"{'arm':<22}{'spec':>8}{'gbF1':>8}{'cgbF1':>8}{'FPR':>8}{'gcount':>8}")
for n in ARMS:
    if n not in rows:
        print(f"{ARMS[n]:<22}{'(pending)':>8}")
        continue
    r = rows[n]
    print(f"{ARMS[n]:<22}{r['spec']:>8.4f}{r['gbF1']:>8.4f}{r['cgbF1']:>8.4f}{r['fpr']:>8.4f}{r['gcount']:>8.3f}")
print(f"\nAnchors: M8 3c gbF1={M8_3C_GBF1} spec={M8_3C_SPEC} | ANNEVO ceiling gbF1={ANNEVO_CEIL}")

print("\n=== M9 CK5 primary_progress_gate ===")
if 0 not in rows:
    print("frozen-L0 missing — STOP.")
    sys.exit(0)
ctrl = rows[0]
print(f"frozen-L0: gbF1={ctrl['gbF1']:.4f} spec={ctrl['spec']:.4f}")
any_pass = False
for n in (2, 4):
    if n not in rows:
        print(f"  L{n}: absent")
        continue
    r = rows[n]
    c1 = r['gbF1'] > ctrl['gbF1']
    c2 = r['gbF1'] > M8_3C_GBF1
    c3 = r['spec'] >= 0.93
    c4 = r['gcount'] >= 0.75
    v = "PASS" if (c1 and c2 and c3 and c4) else "FAIL"
    any_pass = any_pass or v == "PASS"
    print(f"  L{n}: gbF1 {r['gbF1']:.4f} >L0?{c1} >0.7392?{c2} spec {r['spec']:.4f}>=.93?{c3} gc {r['gcount']:.3f}>=.75?{c4} -> {v}")
print(f"\nM9 VERDICT: {'PROGRESS (unfreeze lifts gbF1 -> supplement seeds)' if any_pass else 'NO-LIFT (pivot)'}")
