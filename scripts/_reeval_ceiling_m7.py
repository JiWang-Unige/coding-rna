"""Re-eval REANCHOR-CEILING-ANNEVO-M7 with predictions filtered to TEST chromosomes only.
Fixes the negative-spec bug: ANNEVO predicted on the full subset genome, but eval restricted to
test-chrom genome -> whole-genome predicted-genic bases counted against test-chrom intergenic
denominator (FPR>1). Filter pred GFF to test seqids first. Cheap re-eval, no re-inference.
Run: python scripts/_reeval_ceiling_m7.py
"""
import json, os, subprocess, sys
sys.path.insert(0, ".")
from src.screen_anchor import data as D

EXP = "REANCHOR-CEILING-ANNEVO-M7"
SUB = "outputs/SCREENREF-tiberius_like-ho-s0/eval_subsets"
mets = []
for sp in ["arabidopsis_thaliana", "gallus_gallus"]:
    test_ids = {json.loads(l)["id"] for l in open(f"data/m1_screen/{sp}/split_test.jsonl")}
    raw = f"outputs/{EXP}/predictions/{sp}.gff"
    filt = f"outputs/{EXP}/predictions/{sp}.testfilt.gff"
    D.write_subset_gff(raw, test_ids, filt)
    mj = f"outputs/{EXP}/metrics/{sp}.metrics.json"
    subprocess.run(["python", "scripts/eval_gene_body_mask.py", "--reference-gtf", f"{SUB}/{sp}/reference.gff3",
                    "--prediction-gtf", filt, "--genome-fasta", f"{SUB}/{sp}/genome.fa", "--output-json", mj,
                    "--experiment-id", f"{EXP}_{sp}", "--profile", "screen", "--span-mode", "cds"], check=True)
    m = json.load(open(mj))
    print("  %s: spec %.4f gbF1 %.4f gcount %.3f" % (
        sp, m.get("intergenic_specificity", -1), m.get("gene_body_F1_unconstrained", -1),
        m.get("predicted_gene_count_ratio_vs_reference", -1)))
    mets.append(mj)
subprocess.run(["python", "scripts/aggregate_gene_body_metrics.py", "--metrics", *mets,
                "--output-json", f"outputs/{EXP}/metrics/metrics.json", "--experiment-id", EXP,
                "--profile", "screen"], check=True)
m = json.load(open(f"outputs/{EXP}/metrics/metrics.json"))
print("CEILING(fixed,test-chroms): spec %.4f macro %.4f gbF1 %.4f gcount %.3f" % (
    m.get("intergenic_specificity", -1), m.get("macro_intergenic_specificity", -1),
    m.get("gene_body_F1_unconstrained", -1), m.get("predicted_gene_count_ratio_vs_reference", -1)))
