#!/usr/bin/env python3
"""Validate M1 FASTA/GFF manifest entries before baseline inference."""

import argparse
import hashlib
import json
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit(f"PyYAML is required in the active environment: {exc}")

SPAN_FEATURES = {"CDS", "exon", "intron", "start_codon", "stop_codon"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta_stats(path):
    seqs = 0
    bases = 0
    current = False
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seqs += 1
                current = True
            elif current:
                bases += len(line)
            else:
                raise ValueError("FASTA sequence line before header")
    return {"seqs": seqs, "bases": bases}


def gff_stats(path):
    rows = 0
    features = {}
    genes = set()
    transcripts = set()
    span_rows = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            rows += 1
            feature = cols[2]
            features[feature] = features.get(feature, 0) + 1
            attrs = cols[8]
            for raw in attrs.replace(";", "; ").split(";"):
                part = raw.strip()
                if part.startswith("gene_id "):
                    genes.add(part.split(" ", 1)[1].strip().strip('"'))
                elif part.startswith("ID=gene-") or part.startswith("ID=gene:"):
                    genes.add(part.split("=", 1)[1])
                elif part.startswith("transcript_id "):
                    transcripts.add(part.split(" ", 1)[1].strip().strip('"'))
                elif part.startswith("ID=rna-") or part.startswith("ID=transcript:"):
                    transcripts.add(part.split("=", 1)[1])
            if feature in SPAN_FEATURES:
                span_rows += 1
    return {
        "rows": rows,
        "features": features,
        "gene_ids_seen": len(genes),
        "transcript_ids_seen": len(transcripts),
        "span_rows": span_rows,
    }


def expected_checksum(entry, key):
    checksum = entry.get("checksum")
    if isinstance(checksum, dict):
        return checksum.get(key)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--data-id", action="append", default=[])
    parser.add_argument("--report-json", default=None)
    args = parser.parse_args()

    manifest = yaml.safe_load(Path(args.manifest).read_text(encoding="utf-8"))
    selected = set(args.data_id)
    report = {"status": "pass", "entries": [], "violations": []}

    def violation(data_id, note):
        report["violations"].append({"data_id": data_id, "note": note})

    for entry in manifest.get("pilot_download_queue", []):
        data_id = entry.get("data_id")
        if selected and data_id not in selected:
            continue
        item = {"data_id": data_id, "checks": []}
        genome = Path(entry.get("genome_fasta", ""))
        ref = Path(entry.get("reference_gff_or_gtf", ""))
        if not genome.exists():
            violation(data_id, f"missing genome FASTA: {genome}")
        else:
            actual = sha256(genome)
            item["genome_fasta_sha256"] = actual
            exp = expected_checksum(entry, "genome_fasta_sha256")
            if exp and exp != actual:
                violation(data_id, f"genome checksum mismatch: expected {exp}, got {actual}")
            try:
                item["fasta"] = fasta_stats(genome)
                if item["fasta"]["seqs"] <= 0 or item["fasta"]["bases"] <= 0:
                    violation(data_id, "FASTA is empty")
            except Exception as exc:
                violation(data_id, f"FASTA parse error: {exc}")
        if not ref.exists():
            violation(data_id, f"missing reference GFF/GTF: {ref}")
        else:
            actual = sha256(ref)
            item["reference_gff_or_gtf_sha256"] = actual
            exp = expected_checksum(entry, "reference_gff_or_gtf_sha256")
            if exp and exp != actual:
                violation(data_id, f"reference checksum mismatch: expected {exp}, got {actual}")
            try:
                item["gff"] = gff_stats(ref)
                if item["gff"]["span_rows"] <= 0:
                    violation(data_id, "reference has no CDS/exon/intron/start/stop span rows")
            except Exception as exc:
                violation(data_id, f"GFF/GTF parse error: {exc}")
        item["anchor_eligible"] = entry.get("anchor_eligible")
        item["source_url"] = entry.get("source_url")
        report["entries"].append(item)

    if report["violations"]:
        report["status"] = "fail"
    if args.report_json:
        out = Path(args.report_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["status"] == "pass" else 3)


if __name__ == "__main__":
    main()
