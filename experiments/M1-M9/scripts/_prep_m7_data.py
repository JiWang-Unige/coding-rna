"""REANCHOR-HELDOUT-M7 CK1 data prep (run on baobab login node; file-IO + numpy only).
1) Subset chicken to assembled chromosomes (NC_*) with length <= 20 Mb (drops 6 macrochromosomes
   >20Mb for screen cost budget + 172 unplaced NW_ scaffolds). Preserve full as *.full.*.
   Arabidopsis kept full (7 seqids ~119Mb). anchor & candidate use the SAME subset -> fair.
2) Per species: build 3-class labels, report base-fraction ratios + UTR-rich quantification
   (merged exon_bases - cds_bases), chromosome-level split counts.
3) Emit per-species split_{train,val,test}.jsonl (id=seqid) for the check_data leakage gate.
Run: python scripts/_prep_m7_data.py
"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.screen_anchor import data as D  # noqa: E402

CHICK = os.path.join(ROOT, "data/m1_screen/gallus_gallus")
ARAB = os.path.join(ROOT, "data/m1_screen/arabidopsis_thaliana")
CHICK_MAX_BP = 20_000_000


def merged_bases(gff_path, feat_keep, seqids):
    """Union (merged-interval) base count for a feature type, restricted to seqids."""
    keep = set(seqids)
    ivs = {}
    op = D._open
    with op(gff_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != feat_keep or f[0] not in keep:
                continue
            s0, e1 = int(f[3]) - 1, int(f[4])
            if e1 > s0:
                ivs.setdefault(f[0], []).append((s0, e1))
    total = 0
    for sid, lst in ivs.items():
        lst.sort()
        ce = -1
        for s, e in lst:
            if s > ce:
                total += e - s
                ce = e
            elif e > ce:
                total += e - ce
                ce = e
    return total


def subset_chicken():
    gpath = os.path.join(CHICK, "genome.fa")
    full_g = os.path.join(CHICK, "genome.full.fa")
    full_r = os.path.join(CHICK, "reference.full.gff3")
    rpath = os.path.join(CHICK, "reference.gff3")
    if os.path.exists(full_g):
        print(">>> chicken already subset (genome.full.fa exists); reloading subset seqids")
        seqs = D.read_fasta(gpath)
        return sorted(seqs.keys())
    print(">>> reading FULL chicken genome (~1Gb)...")
    seqs = D.read_fasta(gpath)
    keep = sorted([k for k, v in seqs.items() if k.startswith("NC_") and len(v) <= CHICK_MAX_BP])
    kept_mb = sum(len(seqs[k]) for k in keep) / 1e6
    print(f">>> chicken subset: {len(keep)} NC_ seqids <= 20Mb, {kept_mb:.1f} Mb "
          f"(dropped {len(seqs)-len(keep)} of {len(seqs)})")
    os.rename(gpath, full_g)
    os.rename(rpath, full_r)
    D.write_subset_fasta(seqs, keep, gpath)
    D.write_subset_gff(full_r, keep, rpath)
    print(f">>> wrote subset genome.fa + reference.gff3; full preserved as *.full.*")
    return keep


def species_report(name, sp_dir):
    seqs = D.read_fasta(os.path.join(sp_dir, "genome.fa"))
    seq_lengths = {sid: len(s) for sid, s in seqs.items()}
    labels = D.build_labels(os.path.join(sp_dir, "reference.gff3"), seq_lengths)
    import numpy as np
    allc = np.concatenate([labels[s] for s in seqs])
    tot = len(allc)
    c0 = int((allc == 0).sum()); c1 = int((allc == 1).sum()); c2 = int((allc == 2).sum())
    gff = os.path.join(sp_dir, "reference.gff3")
    exon_b = merged_bases(gff, "exon", seqs.keys())
    cds_b = merged_bases(gff, "CDS", seqs.keys())
    utr_b = max(0, exon_b - cds_b)
    splits = D.assign_splits(list(seqs.keys()))
    from collections import Counter
    sc = Counter(splits.values())
    rep = {
        "species": name, "seqids": len(seqs), "total_bp": tot,
        "class_frac": {"intergenic": round(c0/tot, 4), "CDS": round(c1/tot, 4),
                       "gene_body_nc": round(c2/tot, 4)},
        "exon_bases": exon_b, "cds_bases": cds_b, "utr_bases_est": utr_b,
        "utr_frac_of_exon": round(utr_b/exon_b, 4) if exon_b else 0.0,
        "split_counts": dict(sc),
        "splits": splits,
    }
    # write check_data jsonl (id=seqid) per split
    for sp_name in ("train", "val", "test"):
        rows = [{"id": sid} for sid, t in splits.items() if t == sp_name]
        with open(os.path.join(sp_dir, f"split_{sp_name}.jsonl"), "w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    return rep


def main():
    subset_chicken()
    out = {"chicken_subset_rule": "NC_* assembled chromosomes with length<=20Mb (drops macrochrom>20Mb + NW_ scaffolds); arabidopsis full",
           "species": {}}
    for name, d in [("arabidopsis_thaliana", ARAB), ("gallus_gallus", CHICK)]:
        rep = species_report(name, d)
        out["species"][name] = rep
        cf = rep["class_frac"]
        print(f"\n=== {name} ===")
        print(f"  seqids={rep['seqids']} total={rep['total_bp']/1e6:.1f}Mb "
              f"split={rep['split_counts']}")
        print(f"  class_frac: intergenic={cf['intergenic']} CDS={cf['CDS']} gene_body_nc={cf['gene_body_nc']}")
        print(f"  UTR-rich: exon={rep['exon_bases']/1e6:.1f}Mb CDS={rep['cds_bases']/1e6:.1f}Mb "
              f"UTR_est={rep['utr_bases_est']/1e6:.1f}Mb ({rep['utr_frac_of_exon']*100:.1f}% of exon)")
    op = os.path.join(ROOT, "data/m1_screen/m7_prep_report.json")
    with open(op, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n>>> wrote {op}")


if __name__ == "__main__":
    main()
