#!/usr/bin/env python3
"""Prepare the M13 close-plant diagnostic species: Arabidopsis lyrata.

Downloads RefSeq GCF_000004255.2 v.1.0 genome+GFF, keeps the eight large primary
scaffolds as the screen subset, writes deterministic seqid splits, and records
checksums/provenance. This is data prep only; it does not train or tune a model.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import sys
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.screen_anchor import data as D  # noqa: E402

SPECIES_KEY = "arabidopsis_lyrata"
SPECIES_NAME = "Arabidopsis lyrata subsp. lyrata"
ACCESSION = "GCF_000004255.2"
ASSEMBLY_NAME = "v.1.0"
ANNOTATION_RELEASE = "NCBI Arabidopsis lyrata subsp. lyrata Annotation Release 101"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0"
PREFIX = "GCF_000004255.2_v.1.0"
OUT_DIR = ROOT / "data" / "m1_screen" / SPECIES_KEY
REPORT = ROOT / "reports" / "M13_CLOSE_PLANT_FREEZE.json"

# Assembly report shows the first eight scaffolds are the large A. lyrata primary
# scaffolds (~19-33 Mb each; ~178 Mb total), matching the expected chromosome count.
KEEP_SEQIDS = [
    "NW_003302555.1",
    "NW_003302554.1",
    "NW_003302553.1",
    "NW_003302552.1",
    "NW_003302551.1",
    "NW_003302550.1",
    "NW_003302549.1",
    "NW_003302548.1",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        print(f"skip existing {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    print(f"download {url}")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.rename(path)


def gunzip(src: Path, dst: Path) -> None:
    if dst.exists() and dst.stat().st_size > 0:
        print(f"skip existing {dst}")
        return
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    print(f"decompress {src} -> {dst}")
    with gzip.open(src, "rb") as inp, tmp.open("wb") as out:
        shutil.copyfileobj(inp, out)
    tmp.rename(dst)


def merged_bases(gff_path: Path, feature: str, seqids: set[str]) -> int:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in seqids or fields[2] != feature:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if end > start:
                intervals.setdefault(fields[0], []).append((start, end))
    total = 0
    for vals in intervals.values():
        vals.sort()
        current_end = -1
        for start, end in vals:
            if start > current_end:
                total += end - start
                current_end = end
            elif end > current_end:
                total += end - current_end
                current_end = end
    return total


def write_splits(seqids: list[str], out_dir: Path) -> dict:
    splits = D.assign_splits(seqids)
    for split_name in ("train", "val", "test"):
        rows = [{"id": sid} for sid, val in splits.items() if val == split_name]
        with (out_dir / f"split_{split_name}.jsonl").open("w") as handle:
            for row in rows:
                handle.write(json.dumps(row) + "\n")
    return {"splits": splits, "split_counts": dict(Counter(splits.values()))}


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw_fna = OUT_DIR / f"{PREFIX}_genomic.fna.gz"
    raw_gff = OUT_DIR / f"{PREFIX}_genomic.gff.gz"
    assembly_report = OUT_DIR / f"{PREFIX}_assembly_report.txt"
    readme = OUT_DIR / f"README_Arabidopsis_lyrata_subsp._lyrata_annotation_release_101"

    download(f"{FTP}/{raw_fna.name}", raw_fna)
    download(f"{FTP}/{raw_gff.name}", raw_gff)
    download(f"{FTP}/{assembly_report.name}", assembly_report)
    download(f"{FTP}/{readme.name}", readme)

    full_genome = OUT_DIR / "genome.full.fa"
    full_gff = OUT_DIR / "reference.full.gff3"
    gunzip(raw_fna, full_genome)
    gunzip(raw_gff, full_gff)

    seqs = D.read_fasta(str(full_genome))
    missing = [sid for sid in KEEP_SEQIDS if sid not in seqs]
    if missing:
        raise SystemExit(f"missing expected keep seqids: {missing}")
    D.write_subset_fasta(seqs, KEEP_SEQIDS, str(OUT_DIR / "genome.fa"))
    D.write_subset_gff(str(full_gff), KEEP_SEQIDS, str(OUT_DIR / "reference.gff3"))

    subset_seqs = D.read_fasta(str(OUT_DIR / "genome.fa"))
    seq_lengths = {sid: len(seq) for sid, seq in subset_seqs.items()}
    split_info = write_splits(sorted(subset_seqs), OUT_DIR)
    labels = D.build_labels(str(OUT_DIR / "reference.gff3"), seq_lengths)
    total_bases = sum(seq_lengths.values())
    cds_bases = int(sum((arr == D.CLASS_CDS).sum() for arr in labels.values()))
    nc_bases = int(sum((arr == D.CLASS_GENEBODY_NC).sum() for arr in labels.values()))
    intergenic_bases = total_bases - cds_bases - nc_bases
    exon_bases = merged_bases(OUT_DIR / "reference.gff3", "exon", set(subset_seqs))
    cds_union_bases = merged_bases(OUT_DIR / "reference.gff3", "CDS", set(subset_seqs))

    report = {
        "species_key": SPECIES_KEY,
        "species_name": SPECIES_NAME,
        "accession": ACCESSION,
        "assembly_name": ASSEMBLY_NAME,
        "annotation_release": ANNOTATION_RELEASE,
        "ftp": FTP,
        "status": "downloaded_prepared_frozen_for_m13_diagnostic",
        "selection_rationale": [
            "closest available RefSeq Arabidopsis-relative candidate among checked Brassicaceae options",
            "directly tests whether Arabidopsis-trained fixed model transfers to a near relative",
            "uses top 8 large primary scaffolds to reduce tiny-scaffold noise while keeping screen cost bounded",
        ],
        "caveats": [
            "assembly_level=Scaffold, not chromosome; use for diagnostic only, not final clean claim",
            "do not tune decode/calibration on this species test labels",
            "Brassica rapa GCF_000309985.2 remains chromosome-level backup if scaffold caveat dominates",
        ],
        "paths": {
            "genome_fasta": str(OUT_DIR / "genome.fa"),
            "reference_gff3": str(OUT_DIR / "reference.gff3"),
            "full_genome_fasta": str(full_genome),
            "full_reference_gff3": str(full_gff),
        },
        "subset_rule": {
            "keep_seqids": KEEP_SEQIDS,
            "total_bp": total_bases,
            "seq_lengths": seq_lengths,
            "reason": "top eight large primary scaffolds from assembly report; expected A. lyrata chromosome count is eight",
        },
        "class_fractions": {
            "intergenic": intergenic_bases / total_bases,
            "CDS": cds_bases / total_bases,
            "gene_body_nc": nc_bases / total_bases,
        },
        "annotation_density": {
            "exon_bases": exon_bases,
            "cds_union_bases": cds_union_bases,
            "utr_bases_est": max(0, exon_bases - cds_union_bases),
            "utr_frac_of_exon": (max(0, exon_bases - cds_union_bases) / exon_bases) if exon_bases else 0.0,
        },
        "split": split_info,
        "checksums": {
            "genome_fasta_sha256": sha256(OUT_DIR / "genome.fa"),
            "reference_gff3_sha256": sha256(OUT_DIR / "reference.gff3"),
            "full_genome_fasta_sha256": sha256(full_genome),
            "full_reference_gff3_sha256": sha256(full_gff),
            "raw_genomic_fna_gz_sha256": sha256(raw_fna),
            "raw_genomic_gff_gz_sha256": sha256(raw_gff),
            "assembly_report_sha256": sha256(assembly_report),
        },
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
