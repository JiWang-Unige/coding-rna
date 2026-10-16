#!/usr/bin/env python3
"""M28 development-only context and historical sampler census. No model execution."""
import argparse
import bisect
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from eval_structure_diagnostic import parse_annotation, primary_transcripts
from src.screen_anchor.data import assign_splits

CONFIG = ROOT / "configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml"
WIDTHS = (6144, 12288, 24576, 49152)

def old_sample(descriptors, fraction, cap, seed=0):
    rng = np.random.default_rng(seed)
    if fraction < 1 and descriptors:
        n = max(1, int(round(len(descriptors) * fraction)))
        indices = sorted(rng.choice(len(descriptors), size=n, replace=False))
        descriptors = [descriptors[i] for i in indices]
    return descriptors[:cap] if cap and cap > 0 else descriptors

def uniform_cap(descriptors, cap, seed=0):
    """Prospective control: sample the capped number from the WHOLE pool."""
    n = min(cap, len(descriptors))
    return [descriptors[i] for i in sorted(np.random.default_rng(seed).choice(len(descriptors), n, replace=False))]

def contained(a, b, width, stride, length):
    # At least one regular full window contains the 0-based half-open CDS span.
    lo = max(0, (b - width + stride - 1) // stride)
    hi = min(a // stride, (length - width) // stride)
    return lo <= hi

def union_length(intervals):
    end = -1
    total = 0
    for a, b in sorted(intervals):
        total += max(0, b - max(a, end))
        end = max(end, b)
    return total

def overlapping_ids(txs):
    result = set()
    for strand in ("+", "-"):
        active = []
        for tx in sorted((t for t in txs if t["strand"] == strand), key=lambda t: t["CDS"][0][0]):
            a, b = tx["CDS"][0][0], tx["CDS"][-1][1]
            active = [(end, ident) for end, ident in active if end > a]
            if active:
                result.add(tx["id"])
                result.update(ident for _, ident in active)
            active.append((b, tx["id"]))
    return result

def metadata(config):
    result = {}
    for relpath in config["data"]["development_species"]:
        path = ROOT / relpath
        name = path.name
        with (path / "genome.fa").open() as f:
            ids = [line[1:].split()[0] for line in f if line.startswith(">")]
        assigned = assign_splits(ids)
        for split in ("train", "val"):
            observed = {json.loads(line)["id"] for line in (path / ("split_" + split + ".jsonl")).read_text().splitlines() if line.strip()}
            if observed != {q for q, s in assigned.items() if s == split}:
                raise ValueError("Saved and source-reconstructed split differ: " + name + "/" + split)
        lengths = {}
        with (path / "reference.gff3").open() as f:
            for line in f:
                if line.startswith("##sequence-region "):
                    _, q, a, b = line.split()
                    lengths[q] = int(b) - int(a) + 1
        selected = [q for q in config["data"]["primary_chromosome_seqids"][name] if assigned[q] in ("train", "val")]
        result[name] = {"path": path, "seqids": selected, "splits": {q: assigned[q] for q in selected},
                        "lengths": {q: lengths[q] for q in selected}}
    return result

def descriptors(meta, split, width):
    return [(name, q, start, strand)
            for name, record in meta.items()
            for q in record["seqids"] if record["splits"][q] == split
            for start in range(0, record["lengths"][q] - width + 1, width)
            for strand in ("+", "-")]

def describe_selection(rows):
    counts = Counter((n, q, strand) for n, q, _, strand in rows)
    return {"total": len(rows), "species": dict(Counter(r[0] for r in rows)),
            "chromosome_strand": [{"species": n, "seqid": q, "strand": s, "windows": c}
                                  for (n, q, s), c in sorted(counts.items())]}

def run(species, out):
    config = yaml.safe_load(CONFIG.read_text())
    meta = metadata(config)
    width = config["model"]["window_bp"]
    selected, controls, selection_report = {}, {}, {}
    for split, frac, cap in (("train", .12, 1536), ("val", 1., 768)):
        pool = descriptors(meta, split, width)
        selected[split] = old_sample(pool, frac, cap)
        controls[split] = uniform_cap(pool, cap)
        selection_report[split] = {"pool": describe_selection(pool),
                                   "historical_reconstructed": describe_selection(selected[split]),
                                   "prospective_uniform_same_cap": describe_selection(controls[split])}
    record = meta[species]
    annotation = parse_annotation(record["path"] / "reference.gff3", record["lengths"], protein_coding_only=True)
    all_annotation = parse_annotation(record["path"] / "reference.gff3", record["lengths"])
    complete = primary_transcripts(annotation)
    all_primary = primary_transcripts(annotation, complete_only=False)
    by_q = defaultdict(list)
    for tx in all_primary:
        by_q[tx["seqid"]].append(tx)
    overlap = {q: overlapping_ids(txs) for q, txs in by_q.items()}
    samples = {}
    for split in ("train", "val"):
        for name, rows in (("old", selected[split]), ("uniform", controls[split])):
            by_key = defaultdict(list)
            for n, q, start, strand in rows:
                if n == species:
                    by_key[q, strand].append(start)
            samples[split, name] = by_key
    ledger = []
    for tx in complete:
        q, strand = tx["seqid"], tx["strand"]
        a, b = tx["CDS"][0][0], tx["CDS"][-1][1]
        split = record["splits"][q]
        row = {"species": species, "split": split, "seqid": q, "transcript": tx["id"], "strand": strand,
               "start": a, "end": b, "span": b-a, "CDS_count": len(tx["CDS"]),
               "CDS_bp": sum(y-x for x,y,_ in tx["CDS"]),
               "same_strand_overlap_primary": int(tx["id"] in overlap[q])}
        for name in ("old", "uniform"):
            starts = samples[split, name].get((q, strand), [])
            i = bisect.bisect_right(starts, a) - 1
            row[name + "_sample_full_chain"] = int(i >= 0 and starts[i] + width >= b)
        for w in WIDTHS:
            row["span_fits_" + str(w)] = int(b-a <= w)
            row["tile_contains_" + str(w)] = int(contained(a,b,w,w,record["lengths"][q]))
            row["half_stride_contains_" + str(w)] = int(contained(a,b,w,w//2,record["lengths"][q]))
        ledger.append(row)
    summaries = {}
    for split in ("train", "val"):
        rows = [r for r in ledger if r["split"] == split]
        qids = [q for q in record["seqids"] if record["splits"][q] == split]
        sums = {key: sum(r[key] for r in rows) for key in rows[0] if key.startswith(("old_", "uniform_", "span_fits_", "tile_contains_", "half_stride_"))} if rows else {}
        sampled = [r for r in selected[split] if r[0] == species]
        spans = defaultdict(list)
        for gene in all_annotation["genes"].values():
            spans[gene["seqid"]].append((gene["start"], gene["end"]))
        negative = sum(not any(a < st+width and b > st for a,b in spans[q]) for _,q,st,_ in sampled)
        summaries[split] = {"seqids": qids, "bases": sum(record["lengths"][q] for q in qids),
                           "complete_primary_chains": len(rows),
                           "all_primary_chains_including_partial": sum(t["seqid"] in qids for t in all_primary),
                           "same_strand_overlapping_complete_chains": sum(r["same_strand_overlap_primary"] for r in rows),
                           "CDS_span_quantiles_bp": dict(zip(("p50","p90","p95","p99","max"), np.quantile([r["span"] for r in rows],[.5,.9,.95,.99,1]).tolist())) if rows else {},
                           "counts": sums, "historical_selected_windows": len(sampled),
                           "historical_windows_no_annotated_gene_overlap_either_strand": negative,
                           "all_gene_span_bp": sum(union_length(spans[q]) for q in qids)}
    report = {"experiment": "M28-METHOD-RESTART-CENSUS", "species": species, "status": "development_diagnostic",
              "basis": "current M25R source/config replay, NOT an original saved per-window training trace",
              "no_model_or_training": True, "setaria_access": False, "heldout_sequence_or_labels_used": False,
              "metadata_access": "FASTA headers and GFF sequence-region metadata; retain train/val annotations only",
              "homology_isolation": "not established by chromosome splits",
              "no_innovation_claim": True, "historical_sampler": selection_report, "splits": summaries}
    out.mkdir(parents=True, exist_ok=True)
    with (out / "summary.json").open("x") as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write("\n")
    with (out / "chain_ledger.tsv").open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(ledger[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(ledger)
    print(json.dumps({"species": species, "splits": summaries}, indent=2))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--species", choices=("arabidopsis_thaliana", "oryza_sativa"), required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    run(args.species, args.output_dir)
