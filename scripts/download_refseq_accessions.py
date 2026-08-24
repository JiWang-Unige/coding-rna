#!/usr/bin/env python3
"""Download RefSeq genome FASTA + GFF3 for accessions in a YAML manifest."""

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit(f"PyYAML is required in the active environment: {exc}")


ASSEMBLY_SUMMARY = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
)


def run(cmd, *, check=True):
    return subprocess.run(cmd, text=True, capture_output=True, check=check)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    cmd = ["curl", "-L", "--fail", "--retry", "3", "--retry-delay", "5", "-o", str(dest), url]
    run(cmd)


def gunzip_to(src_gz, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src_gz, "rb") as src, open(dest, "wb") as out:
        shutil.copyfileobj(src, out)


def load_assembly_summary(cache_path):
    if not cache_path.exists() or cache_path.stat().st_size == 0:
        download(ASSEMBLY_SUMMARY, cache_path)
    rows = {}
    with open(cache_path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 20:
                continue
            accession = cols[0]
            rows[accession] = {
                "accession": accession,
                "organism_name": cols[7],
                "assembly_name": cols[15],
                "ftp_path": cols[19],
            }
    return rows


def file_prefix(ftp_path):
    return ftp_path.rstrip("/").split("/")[-1]


def normalize_entry(entry, assembly_rows, cache_dir):
    accession = entry.get("accession")
    if not accession or not accession.startswith("GCF_"):
        return None
    if accession not in assembly_rows:
        raise RuntimeError(f"{accession} not found in RefSeq assembly summary")
    row = assembly_rows[accession]
    ftp = row["ftp_path"]
    prefix = file_prefix(ftp)
    raw_dir = cache_dir / accession
    genomic_gz = raw_dir / f"{prefix}_genomic.fna.gz"
    gff_gz = raw_dir / f"{prefix}_genomic.gff.gz"
    download(f"{ftp}/{prefix}_genomic.fna.gz", genomic_gz)
    download(f"{ftp}/{prefix}_genomic.gff.gz", gff_gz)

    genome_path = Path(entry["genome_fasta"])
    ref_path = Path(entry["reference_gff_or_gtf"])
    gunzip_to(genomic_gz, genome_path)
    gunzip_to(gff_gz, ref_path)

    checks = {
        "genome_fasta_sha256": sha256(genome_path),
        "reference_gff_or_gtf_sha256": sha256(ref_path),
        "raw_genomic_fna_gz_sha256": sha256(genomic_gz),
        "raw_genomic_gff_gz_sha256": sha256(gff_gz),
    }
    entry["source_url"] = ftp
    entry["assembly_name"] = row["assembly_name"]
    entry["organism_name"] = row["organism_name"]
    entry["checksum"] = checks
    entry["download_status"] = "downloaded"
    return {
        "data_id": entry.get("data_id"),
        "accession": accession,
        "ftp_path": ftp,
        "genome_fasta": str(genome_path),
        "reference_gff_or_gtf": str(ref_path),
        **checks,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--data-id", action="append", default=[])
    parser.add_argument("--cache-dir", default="data/_refseq_cache")
    parser.add_argument("--report-json", default=None)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    summary = load_assembly_summary(cache_dir / "assembly_summary_refseq.txt")

    selected = set(args.data_id)
    reports = []
    for entry in manifest.get("pilot_download_queue", []):
        if selected and entry.get("data_id") not in selected:
            continue
        rep = normalize_entry(entry, summary, cache_dir)
        if rep:
            reports.append(rep)

    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    if args.report_json:
        out = Path(args.report_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(reports, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
