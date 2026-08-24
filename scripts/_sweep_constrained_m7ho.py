"""REANCHOR-HELDOUT-M7 CK6 — VAL-chosen constrained-param sweep for the held-out candidate
(no test leakage), to land gene_count in the two-sided band [1.0,1.25] on the UTR-rich cross-clade
species. Torch-free (copies constrained_decode; labels_to_cds_gff + eval are numpy/stdlib).
Loads RAW (pre-constrained) per-seqid preds from outputs/<run>/raw_pred/{val,test}_{sp}.npz.
Held-out version of scripts/_sweep_constrained_m6.py (SP + runs changed).
Run: python scripts/_sweep_constrained_m7ho.py
"""
import argparse, json, os, subprocess, sys, tempfile
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.screen_anchor.gff_io import labels_to_cds_gff      # numpy-only  # noqa: E402

CLASS_CDS, CLASS_GENEBODY_NC = 1, 2
SP = ["arabidopsis_thaliana", "gallus_gallus"]
PY = sys.executable


def constrained_decode(pred_by_seqid, min_cds_len=30, max_fill_gap=20):
    out = {}
    for sid, arr in pred_by_seqid.items():
        a = arr.copy(); n = len(a)
        is_cds = a == CLASS_CDS
        i = 0
        while i < n:
            if is_cds[i]:
                j = i
                while j < n and is_cds[j]:
                    j += 1
                if (j - i) < min_cds_len:
                    a[i:j] = CLASS_GENEBODY_NC
                i = j
            else:
                i += 1
        gene = a > 0
        i = 0
        while i < n:
            if not gene[i]:
                j = i
                while j < n and not gene[j]:
                    j += 1
                if 0 < i and j < n and (j - i) <= max_fill_gap:
                    a[i:j] = CLASS_GENEBODY_NC
                i = j
            else:
                i += 1
        out[sid] = a
    return out


def decode_eval(run, split, mfg, mcl, tmp):
    subset = "val_eval_subsets" if split == "val" else "eval_subsets"
    sj = []
    for sp in SP:
        npz = f"{ROOT}/outputs/{run}/raw_pred/{split}_{sp}.npz"
        ref = f"{ROOT}/outputs/{run}/{subset}/{sp}/reference.gff3"
        gen = f"{ROOT}/outputs/{run}/{subset}/{sp}/genome.fa"
        if not all(os.path.exists(x) for x in (npz, ref, gen)):
            return None
        d = np.load(npz, allow_pickle=True)
        cd = constrained_decode({k: d[k] for k in d.files}, min_cds_len=mcl, max_fill_gap=mfg)
        gff = os.path.join(tmp, f"{run}_{split}_{sp}.gff"); labels_to_cds_gff(cd, gff)
        oj = os.path.join(tmp, f"{run}_{split}_{sp}.json")
        subprocess.run([PY, f"{ROOT}/scripts/eval_gene_body_mask.py", "--reference-gtf", ref,
                        "--prediction-gtf", gff, "--genome-fasta", gen, "--output-json", oj,
                        "--experiment-id", f"{run}_{split}", "--profile", "screen", "--span-mode", "cds"],
                       check=True)
        sj.append(oj)
    agg = os.path.join(tmp, f"{run}_{split}_AGG.json")
    subprocess.run([PY, f"{ROOT}/scripts/aggregate_gene_body_metrics.py", "--metrics", *sj,
                    "--output-json", agg, "--experiment-id", f"{run}_{split}", "--profile", "screen"], check=True)
    return json.load(open(agg))


def mean_over_runs(runs, split, mfg, mcl, tmp, key):
    vals = []
    for r in runs:
        a = decode_eval(r, split, mfg, mcl, tmp)
        if a:
            vals.append(a[key])
    return float(np.mean(vals)) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", default=[f"FP-FRAGFIX-CONSTR-ho-s{s}" for s in range(5)])
    ap.add_argument("--mfg", nargs="+", type=int, default=[20, 40, 60, 100, 150])
    ap.add_argument("--mcl", nargs="+", type=int, default=[10, 20, 30, 60, 90])
    ap.add_argument("--anchor", type=float, default=0.8054, help="held-out 3-seed anchor spec (strict-exceed bar)")
    ap.add_argument("--lower-band", type=float, default=1.0,
                    help="TWO-SIDED gene_count band: eligible iff lower-band <= val_gcount <= 1.25, then max val_spec.")
    args = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="m7hosweep_")
    print(f"runs={args.runs}\ngrid mfg={args.mfg} mcl={args.mcl}\ntmp={tmp}\n")
    print("=== VAL grid (choose params here; NO test) ===")
    print(f"{'mfg':>5}{'mcl':>5}{'val_spec':>10}{'val_gcount':>12}{'in_band':>9}")
    rows = []
    for mfg in args.mfg:
        for mcl in args.mcl:
            vs = mean_over_runs(args.runs, "val", mfg, mcl, tmp, "intergenic_specificity")
            vg = mean_over_runs(args.runs, "val", mfg, mcl, tmp, "predicted_gene_count_ratio_vs_reference")
            if vs is None:
                print(f"{mfg:>5}{mcl:>5}  MISSING raw_pred"); return
            inb = args.lower_band <= vg <= 1.25
            rows.append((mfg, mcl, vs, vg, inb))
            print(f"{mfg:>5}{mcl:>5}{vs:>10.4f}{vg:>12.3f}{('Y' if inb else 'n'):>9}")
    elig = [r for r in rows if r[4]]
    if not elig:
        # fallback: closest to 1.0 from the grid (report, do not silently fail)
        best = min(rows, key=lambda r: abs(r[3] - 1.0))
        print(f"\nNo grid point in band [{args.lower_band},1.25]; FALLBACK closest-to-1.0: mfg={best[0]} mcl={best[1]} val_gcount={best[3]:.3f}")
    else:
        best = max(elig, key=lambda r: r[2])
        print(f"\n(two-sided band [{args.lower_band},1.25]; {len(elig)} eligible)")
    mfg, mcl = best[0], best[1]
    print(f"CHOSEN (on VAL): max_fill_gap={mfg} min_cds_len={mcl} (val_spec {best[2]:.4f}, val_gcount {best[3]:.3f})")
    print("=== APPLY to TEST (5-seed, mean+-std) ===")
    for k, lab in [("intergenic_specificity", "spec"), ("macro_intergenic_specificity", "macro"),
                   ("gene_body_F1_unconstrained", "gbF1"), ("predicted_gene_count_ratio_vs_reference", "gcount")]:
        per = [decode_eval(r, "test", mfg, mcl, tmp) for r in args.runs]
        per = [a[k] for a in per if a]
        print(f"  test {lab:7s} mean={np.mean(per):.4f} +- {np.std(per):.4f} per-seed={[round(x,4) for x in per]}")
    print(f"\nGATE: test spec > held-out anchor {args.anchor} AND gbF1>=0.5276 AND macro>=0.7804 AND gcount in [1.0,1.25]")


if __name__ == "__main__":
    main()
