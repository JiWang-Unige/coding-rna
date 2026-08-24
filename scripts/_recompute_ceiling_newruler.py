"""Recompute pretrained_ceiling under the NEW full-transcript intergenic ruler (R6).
Re-evaluate the whole-genome pretrained Tiberius/Helixer/ANNEVO prediction GFFs on the SAME
held-out test subset as the anchor (filter predictions to test seqids), --span-mode cds, then
aggregate base-weighted + macro intergenic_specificity. ceiling = max over the 3 tools.
Pure stdlib eval (no GPU). Run: python3 scripts/_recompute_ceiling_newruler.py
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
SPECIES = ["saccharomyces_cerevisiae", "drosophila_melanogaster"]
SUBSET = ROOT / "outputs/SCREENREF-tiberius_like-s0/eval_subsets"   # test-subset ref+genome (anchor split)
W = Path("/tmp/ceiling_newruler"); W.mkdir(parents=True, exist_ok=True)

TOOLS = {
    "Tiberius": {sp: ROOT / f"outputs/BASE-TIBERIUS-PILOT-M1/predictions/{sp}.gtf" for sp in SPECIES},
    "Helixer":  {sp: ROOT / f"outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/predictions/{sp}.gff3" for sp in SPECIES},
    "ANNEVO":   {sp: ROOT / f"outputs/BASE-ANNEVO-SAC-DMEL-SMOKE-M1/predictions/{sp}.gff" for sp in SPECIES},
}


def seqids_of(fasta):
    out = set()
    for ln in open(fasta):
        if ln.startswith(">"):
            out.add(ln[1:].split()[0])
    return out


def main():
    print("tool       | spec_bw | spec_macro | FPR_bw | gbF1_bw")
    print("-" * 60)
    ceiling = {}
    for tool, preds in TOOLS.items():
        sp_jsons = []
        ok = True
        for sp in SPECIES:
            ref = SUBSET / sp / "reference.gff3"
            genome = SUBSET / sp / "genome.fa"
            pred = preds[sp]
            if not pred.exists():
                print(f"{tool}: missing {pred}"); ok = False; break
            keep = seqids_of(genome)
            fp = W / f"{tool}_{sp}.gff"
            with open(pred) as fi, open(fp, "w") as fo:
                for ln in fi:
                    if ln.startswith("#") or ln.split("\t", 1)[0] in keep:
                        fo.write(ln)
            outj = W / f"{tool}_{sp}.json"
            subprocess.run([PY, str(ROOT / "scripts/eval_gene_body_mask.py"),
                            "--reference-gtf", str(ref), "--prediction-gtf", str(fp),
                            "--genome-fasta", str(genome), "--output-json", str(outj),
                            "--experiment-id", f"CEIL_{tool}_{sp}", "--profile", "screen",
                            "--span-mode", "cds"], check=True)
            sp_jsons.append(str(outj))
        if not ok:
            continue
        agg = W / f"{tool}_AGG.json"
        subprocess.run([PY, str(ROOT / "scripts/aggregate_gene_body_metrics.py"),
                        "--metrics", *sp_jsons, "--output-json", str(agg),
                        "--experiment-id", f"CEIL_{tool}", "--profile", "screen"], check=True)
        d = json.load(open(agg))
        ceiling[tool] = d
        print("%-10s | %.4f | %.4f | %.4f | %.4f" % (
            tool, d["intergenic_specificity"], d.get("macro_intergenic_specificity") or 0,
            d["intergenic_FPR"], d["gene_body_F1_unconstrained"]))
    print("-" * 60)
    if ceiling:
        best = max(ceiling, key=lambda t: ceiling[t]["intergenic_specificity"])
        b = ceiling[best]
        print("CEILING (max intergenic_specificity) = %s  spec_bw=%.4f macro=%.4f gbF1=%.4f" % (
            best, b["intergenic_specificity"], b.get("macro_intergenic_specificity") or 0,
            b["gene_body_F1_unconstrained"]))
    json.dump(ceiling, open(W / "ceiling_summary.json", "w"), indent=2, default=float)


if __name__ == "__main__":
    main()
