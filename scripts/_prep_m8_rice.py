"""M8-CK1 (clean held-out): subset rice (Oryza sativa, 386Mb full) to ~150Mb for screen cost,
mirroring the chicken rule. Reproducible: assembled NC_ chromosomes, ascending by size, cumulative
<=150Mb. Preserve full as *.full.*. Then emit split jsonl + run nothing else (MC labels already
verified by _verify_mc_labels.py). Run: python scripts/_prep_m8_rice.py
"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.screen_anchor import data as D

RICE = os.path.join(ROOT, "data/m1_screen/oryza_sativa")
BUDGET_MB = 150.0


def subset_rice():
    g = os.path.join(RICE, "genome.fa"); full_g = os.path.join(RICE, "genome.full.fa")
    r = os.path.join(RICE, "reference.gff3"); full_r = os.path.join(RICE, "reference.full.gff3")
    if os.path.exists(full_g):
        print(">>> rice already subset; reloading"); return sorted(D.read_fasta(g).keys())
    seqs = D.read_fasta(g)
    nc = sorted([(k, len(v)) for k, v in seqs.items() if k.startswith("NC_")], key=lambda x: x[1])
    keep, cum = [], 0.0
    for k, L in nc:
        if cum + L / 1e6 > BUDGET_MB and keep:
            break
        keep.append(k); cum += L / 1e6
    keep = sorted(keep)
    print(f">>> rice subset: {len(keep)} NC_ chromosomes, {cum:.1f}Mb (of {len(seqs)} seqids, "
          f"{sum(len(v) for v in seqs.values())/1e6:.1f}Mb full)")
    os.rename(g, full_g); os.rename(r, full_r)
    D.write_subset_fasta(seqs, keep, g)
    D.write_subset_gff(full_r, keep, r)
    return keep


def main():
    subset_rice()
    seqs = D.read_fasta(os.path.join(RICE, "genome.fa"))
    splits = D.assign_splits(list(seqs.keys()))
    from collections import Counter
    print(">>> rice split:", dict(Counter(splits.values())))
    for sp_name in ("train", "val", "test"):
        rows = [{"id": sid} for sid, t in splits.items() if t == sp_name]
        with open(os.path.join(RICE, f"split_{sp_name}.jsonl"), "w") as fh:
            for rrow in rows:
                fh.write(json.dumps(rrow) + "\n")
    print(">>> wrote split_{train,val,test}.jsonl for oryza_sativa")


if __name__ == "__main__":
    main()
