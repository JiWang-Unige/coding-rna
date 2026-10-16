import numpy as np, os
d = "outputs/FP-FRAGFIX-RP-SMOKE/raw_pred"
print("raw_pred files:", sorted(os.listdir(d)) if os.path.isdir(d) else "MISSING")
for f in ["val_saccharomyces_cerevisiae.npz", "test_saccharomyces_cerevisiae.npz",
          "val_drosophila_melanogaster.npz", "test_drosophila_melanogaster.npz"]:
    p = os.path.join(d, f)
    if os.path.exists(p):
        z = np.load(p, allow_pickle=True); ks = z.files; a = z[ks[0]]
        print(f"  {f}: {len(ks)} seqids; {ks[0]} shape={a.shape} classes={sorted(set(a.tolist()))[:4]}")
    else:
        print(f"  {f}: MISSING")
vs = "outputs/FP-FRAGFIX-RP-SMOKE/val_eval_subsets"
print("val_eval_subsets:", os.listdir(vs) if os.path.isdir(vs) else "MISSING")
