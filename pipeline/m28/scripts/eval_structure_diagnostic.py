#!/usr/bin/env python3
"""M24: same-scope structural diagnostics for saved direct-annotation artifacts."""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import eval_gene_body_mask as gene_body


SPECIES = ("arabidopsis_thaliana", "oryza_sativa")
METHODS = {
    "M19_s0": ("outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0/predictions", ".gff", "candidate"),
    "M19_s1": ("outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/predictions", ".gff", "candidate"),
    "M8_3C_s0": ("outputs/M8-3C-CAND-s0/predictions", ".gff", "candidate"),
    "M8_3C_s2": ("outputs/M8-3C-CAND-s2/predictions", ".gff", "candidate"),
    "M8_3C_s4": ("outputs/M8-3C-CAND-s4/predictions", ".gff", "candidate"),
    "ANNEVO": ("outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/predictions", ".gff", "baseline_coding"),
    "Helixer": ("outputs/M12B-SAMEPANEL-BASELINES-HELIXER/predictions", ".gff3", "baseline_full"),
    "Tiberius": ("outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/predictions", ".gtf", "baseline_coding"),
}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "rRNA", "tRNA"}
NA = "not_applicable"


def parse_attrs(text):
    attrs = {}
    for item in text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item.split(None, 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def is_partial(attrs):
    return attrs.get("partial", "").lower() == "true" or "start_range" in attrs or "end_range" in attrs


def parse_annotation(path, lengths, protein_coding_only=False):
    records = []
    retained = 0
    out_of_bounds = 0
    for line_no, line in enumerate(open(path), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"{path}:{line_no}: expected 9 columns")
        seqid, _source, feature, start_s, end_s, _score, strand, phase, attrs_s = fields
        if seqid not in lengths:
            continue
        start, end = int(start_s) - 1, int(end_s)
        retained += 1
        if start < 0 or end <= start or end > lengths[seqid]:
            out_of_bounds += 1
            continue
        records.append((seqid, feature, start, end, strand, phase, parse_attrs(attrs_s), line_no))
    if out_of_bounds:
        raise ValueError(f"{path}: {out_of_bounds} out-of-bounds records")

    genes = {}
    transcripts = {}
    for seqid, feature, start, end, strand, _phase, attrs, _line_no in records:
        if feature == "gene":
            gene_id = attrs.get("ID") or attrs.get("gene_id")
            if not gene_id:
                raise ValueError(f"{path}: gene without ID")
            biotype = attrs.get("gene_biotype") or attrs.get("biotype") or attrs.get("gene_type")
            genes[gene_id] = {
                "id": gene_id, "seqid": seqid, "strand": strand, "start": start, "end": end,
                "coding": biotype == "protein_coding" if protein_coding_only else True,
                "partial": is_partial(attrs), "transcripts": [],
            }
        elif feature in TRANSCRIPT_TYPES:
            tx_id = attrs.get("ID") or attrs.get("transcript_id")
            gene_id = attrs.get("gene_id") or attrs.get("Parent")
            if not tx_id or not gene_id or "," in gene_id:
                raise ValueError(f"{path}: ambiguous transcript parent")
            transcripts[tx_id] = {
                "id": tx_id, "gene_id": gene_id, "seqid": seqid, "strand": strand,
                "partial": is_partial(attrs), "CDS": [], "exon": [],
            }

    for seqid, feature, start, end, strand, phase, attrs, line_no in records:
        if feature not in {"CDS", "exon", "intron"}:
            continue
        parents = []
        if attrs.get("Parent"):
            parents = attrs["Parent"].split(",")
        elif attrs.get("transcript_id"):
            parents = [attrs["transcript_id"]]
        elif attrs.get("gene_id"):
            parents = [attrs["gene_id"]]
        if len(parents) != 1:
            raise ValueError(f"{path}:{line_no}: ambiguous feature parent")
        tx_id = parents[0]
        if tx_id not in transcripts:
            gene_id = attrs.get("gene_id") or tx_id
            transcripts[tx_id] = {
                "id": tx_id, "gene_id": gene_id, "seqid": seqid, "strand": strand,
                "partial": is_partial(attrs), "CDS": [], "exon": [],
            }
        tx = transcripts[tx_id]
        if tx["seqid"] != seqid or tx["strand"] != strand:
            raise ValueError(f"{path}:{line_no}: inconsistent transcript coordinates")
        if feature in {"CDS", "exon"}:
            tx[feature].append((start, end, phase if feature == "CDS" else "."))

    for tx in transcripts.values():
        gene_id = tx["gene_id"]
        if gene_id not in genes:
            if protein_coding_only:
                continue
            spans = tx["CDS"] or tx["exon"]
            if not spans:
                continue
            genes[gene_id] = {
                "id": gene_id, "seqid": tx["seqid"], "strand": tx["strand"],
                "start": min(x[0] for x in spans), "end": max(x[1] for x in spans),
                "coding": True, "partial": False, "transcripts": [],
            }
        genes[gene_id]["transcripts"].append(tx["id"])

    if protein_coding_only:
        keep_genes = {gene_id for gene_id, gene in genes.items() if gene["coding"]}
        genes = {gene_id: gene for gene_id, gene in genes.items() if gene_id in keep_genes}
        transcripts = {tx_id: tx for tx_id, tx in transcripts.items() if tx["gene_id"] in keep_genes}
    for tx in transcripts.values():
        tx["CDS"].sort()
        tx["exon"].sort()
    return {"genes": genes, "transcripts": transcripts, "retained_records": retained, "out_of_bounds": 0}


def primary_transcripts(annotation, complete_only=True):
    selected = []
    for gene in annotation["genes"].values():
        if complete_only and gene["partial"]:
            continue
        choices = []
        for tx_id in gene["transcripts"]:
            tx = annotation["transcripts"].get(tx_id)
            if tx and tx["CDS"] and (not complete_only or not tx["partial"]):
                choices.append(tx)
        if choices:
            selected.append(min(choices, key=lambda tx: (-sum(end - start for start, end, _ in tx["CDS"]), tx["id"])))
    return selected


def intervals(transcripts, feature="CDS", strand=True):
    values = set()
    for tx in transcripts:
        for start, end, phase in tx[feature]:
            key = (tx["seqid"], start, end, tx["strand"]) if strand else (tx["seqid"], start, end)
            values.add(key + ((phase,) if feature == "CDS" and strand else ()))
    return values


def chains(transcripts, strand=True, feature="CDS"):
    values = set()
    for tx in transcripts:
        parts = tuple((start, end) for start, end, _ in tx[feature])
        if parts:
            values.add((tx["seqid"], tx["strand"], parts) if strand else (tx["seqid"], parts))
    return values


def introns(transcript, feature="CDS"):
    parts = transcript[feature]
    return tuple((parts[i][1], parts[i + 1][0]) for i in range(len(parts) - 1) if parts[i][1] < parts[i + 1][0])


def intron_chains(transcripts, strand=True, feature="CDS"):
    values = set()
    for tx in transcripts:
        chain = introns(tx, feature)
        if chain:
            values.add((tx["seqid"], tx["strand"], chain) if strand else (tx["seqid"], chain))
    return values


def intron_intervals(transcripts, strand=True, feature="CDS"):
    values = set()
    for tx in transcripts:
        for start, end in introns(tx, feature):
            values.add((tx["seqid"], start, end, tx["strand"]) if strand
                       else (tx["seqid"], start, end))
    return values


def prf(predicted, reference):
    matched = len(predicted & reference)
    precision = matched / len(predicted) if predicted else 0.0
    recall = matched / len(reference) if reference else 0.0
    return {"matched": matched, "predicted": len(predicted), "reference": len(reference),
            "precision": precision, "recall": recall,
            "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0}


def boundary_metric(predicted_intervals, reference_intervals, side, tolerance):
    pred = defaultdict(list)
    ref = defaultdict(list)
    index = 1 if side == "left" else 2
    for item in predicted_intervals:
        pred[item[0]].append(item[index])
    for item in reference_intervals:
        ref[item[0]].append(item[index])
    matched = 0
    for seqid in set(pred) | set(ref):
        a, b = sorted(pred[seqid]), sorted(ref[seqid])
        i = j = 0
        while i < len(a) and j < len(b):
            if abs(a[i] - b[j]) <= tolerance:
                matched += 1
                i += 1
                j += 1
            elif a[i] < b[j]:
                i += 1
            else:
                j += 1
    precision = matched / len(predicted_intervals) if predicted_intervals else 0.0
    recall = matched / len(reference_intervals) if reference_intervals else 0.0
    return {"matched": matched, "precision": precision, "recall": recall}


def overlap_degrees(a_transcripts, b_transcripts):
    events = []
    for side, transcripts in (("a", a_transcripts), ("b", b_transcripts)):
        for tx in transcripts:
            parts = tx["CDS"]
            if parts:
                events.append((tx["seqid"], parts[0][0], 1, side, tx["id"]))
                events.append((tx["seqid"], parts[-1][1], 0, side, tx["id"]))
    active = {"a": set(), "b": set()}
    degree = {"a": Counter(), "b": Counter()}
    last_seqid = None
    for seqid, _position, kind, side, item_id in sorted(events):
        if seqid != last_seqid:
            active = {"a": set(), "b": set()}
            last_seqid = seqid
        other = "b" if side == "a" else "a"
        if kind == 0:
            active[side].remove(item_id)
        else:
            for other_id in active[other]:
                degree[side][item_id] += 1
                degree[other][other_id] += 1
            active[side].add(item_id)
    return {
        "predicted_groups_overlapping_multiple_reference_transcripts": sum(v > 1 for v in degree["a"].values()),
        "reference_transcripts_overlapping_multiple_predicted_groups": sum(v > 1 for v in degree["b"].values()),
    }


def transcript_spans(annotation, full=False):
    by_seqid = defaultdict(list)
    for tx in annotation["transcripts"].values():
        parts = tx["exon"] if full and tx["exon"] else tx["CDS"]
        if parts:
            by_seqid[tx["seqid"]].append((parts[0][0], parts[-1][1]))
    return {seqid: gene_body.merge_intervals(values) for seqid, values in by_seqid.items()}


def gene_body_metrics(reference, prediction, lengths):
    ref_cds = transcript_spans(reference)
    ref_full = transcript_spans(reference, full=True)
    pred_cds = transcript_spans(prediction)
    ref_len = gene_body.interval_length(ref_cds)
    pred_len = gene_body.interval_length(pred_cds)
    overlap = gene_body.intersection_length(ref_cds, pred_cds)
    precision = overlap / pred_len if pred_len else 0.0
    recall = overlap / ref_len if ref_len else 0.0
    ref_full_len = gene_body.interval_length(ref_full)
    pred_in_full = gene_body.intersection_length(pred_cds, ref_full)
    intergenic = sum(lengths.values()) - ref_full_len
    ref_gene_count = sum(any(reference["transcripts"].get(tx_id, {}).get("CDS") for tx_id in gene["transcripts"])
                         for gene in reference["genes"].values())
    pred_gene_count = sum(any(prediction["transcripts"].get(tx_id, {}).get("CDS") for tx_id in gene["transcripts"])
                          for gene in prediction["genes"].values())
    return {
        "gene_body_F1": gene_body.f1(precision, recall),
        "intergenic_FPR": (pred_len - pred_in_full) / intergenic if intergenic else 0.0,
        "predicted_gene_count_ratio": pred_gene_count / ref_gene_count if ref_gene_count else 0.0,
    }


def exact_structure_metrics(reference, prediction, capability):
    ref_tx = primary_transcripts(reference)
    pred_tx = primary_transcripts(prediction)
    ref_coord = intervals(ref_tx, strand=False)
    pred_coord = intervals(pred_tx, strand=False)
    result = {
        "exact_CDS_run": prf(pred_coord, ref_coord),
        "boundaries": {},
        "pseudo_CDS_chain": prf(chains(pred_tx, strand=False), chains(ref_tx, strand=False)),
        "overlap_fragmentation": overlap_degrees(pred_tx, ref_tx),
    }
    for tolerance in (0, 1, 3, 6):
        result["boundaries"][str(tolerance)] = {
            side: boundary_metric(pred_coord, ref_coord, side, tolerance) for side in ("left", "right")
        }
    if capability == "candidate":
        result.update({
            "strand_aware_CDS": NA, "strand_accuracy": NA, "phase_accuracy": NA,
            "exact_CDS_chain": NA, "exact_intron": NA, "exact_intron_chain": NA,
            "exact_coding_transcript": NA, "exact_coding_gene": NA,
            "exact_full_exon_transcript": NA, "exact_full_exon_gene": NA,
        })
        return result

    ref_strand = {(a, b, c, d) for a, b, c, d, _phase in intervals(ref_tx)}
    pred_strand = {(a, b, c, d) for a, b, c, d, _phase in intervals(pred_tx)}
    ref_phase = {(a, b, c, d): phase for a, b, c, d, phase in intervals(ref_tx)}
    pred_phase = {(a, b, c, d): phase for a, b, c, d, phase in intervals(pred_tx)}
    comparable = set(ref_phase) & set(pred_phase)
    represented = [key for key in comparable if ref_phase[key] != "." and pred_phase[key] != "."]
    coordinate_matches = set((a, b, c) for a, b, c, _strand in ref_strand) & set((a, b, c) for a, b, c, _strand in pred_strand)
    ref_strands = {(a, b, c): strand for a, b, c, strand in ref_strand}
    pred_strands = {(a, b, c): strand for a, b, c, strand in pred_strand}
    result.update({
        "strand_aware_CDS": prf(pred_strand, ref_strand),
        "strand_accuracy": (sum(pred_strands[key] == ref_strands[key] for key in coordinate_matches) /
                            len(coordinate_matches) if coordinate_matches else NA),
        "phase_accuracy": (sum(pred_phase[key] == ref_phase[key] for key in represented) / len(represented)
                           if represented else NA),
        "exact_CDS_chain": prf(chains(pred_tx), chains(ref_tx)),
        "exact_intron": prf(intron_intervals(pred_tx), intron_intervals(ref_tx)),
        "exact_intron_chain": prf(intron_chains(pred_tx), intron_chains(ref_tx)),
        "exact_coding_transcript": prf(chains(pred_tx), chains(ref_tx)),
        "exact_coding_gene": prf(chains(pred_tx), chains(ref_tx)),
    })
    if capability == "baseline_full":
        result["exact_full_exon_transcript"] = prf(chains(pred_tx, feature="exon"), chains(ref_tx, feature="exon"))
        result["exact_full_exon_gene"] = result["exact_full_exon_transcript"]
    else:
        result["exact_full_exon_transcript"] = NA
        result["exact_full_exon_gene"] = NA
    return result


def truth_mask(length, transcripts, element):
    mask = np.zeros(length, dtype=bool)
    for tx in transcripts:
        if element == "exon":
            regions = [(start, end) for start, end, _ in tx["exon"]]
        elif element == "intron":
            regions = list(introns(tx, "exon" if tx["exon"] else "CDS"))
        else:
            regions = []
            for start, end in introns(tx, "exon" if tx["exon"] else "CDS"):
                if element == "splice_donor":
                    regions.append((start, min(start + 2, end)) if tx["strand"] == "+" else (max(start, end - 2), end))
                else:
                    regions.append((max(start, end - 2), end) if tx["strand"] == "+" else (start, min(start + 2, end)))
        for start, end in regions:
            mask[start:end] = True
    return mask


def score_histogram(scores, truth):
    bits = np.asarray(scores, dtype=np.float16).view(np.uint16)
    positives = np.bincount(bits[truth], minlength=65536)
    negatives = np.bincount(bits[~truth], minlength=65536)
    tp = np.cumsum(positives[::-1])
    fp = np.cumsum(negatives[::-1])
    total_positive = int(truth.sum())
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp, dtype=float), where=(tp + fp) != 0)
    recall = tp / total_positive if total_positive else np.zeros_like(tp, dtype=float)
    ap = float(np.sum(np.diff(np.r_[0.0, recall]) * precision))
    f1 = np.divide(2 * precision * recall, precision + recall,
                   out=np.zeros_like(precision), where=(precision + recall) != 0)
    best = int(np.argmax(f1))
    threshold_bits = np.uint16(65535 - best)
    threshold = float(np.array([threshold_bits], dtype=np.uint16).view(np.float16)[0])
    return ap, threshold, float(f1[best])


def score_at_threshold(scores, truth, threshold):
    predicted = scores >= threshold
    tp = int(np.count_nonzero(predicted & truth))
    precision = tp / int(predicted.sum()) if predicted.any() else 0.0
    recall = tp / int(truth.sum()) if truth.any() else 0.0
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def segmentnt_metrics(root, species):
    cache_path = root / "outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species" / f"{species}.npz"
    data = np.load(cache_path, allow_pickle=False)
    feature_names = data["features"].astype(str).tolist()
    splits = dict(zip(data["seqids"].astype(str), data["splits"].astype(str)))
    val_seqid = next(seqid for seqid, split in splits.items() if split == "val")
    test_seqid = next(seqid for seqid, split in splits.items() if split == "test")
    reference_path = root / "data/m1_screen" / species / "reference.gff3"
    lengths = gene_body.fasta_lengths(root / "data/m1_screen" / species / "genome.fa")
    scored_lengths = {seqid: lengths[seqid] for seqid in (val_seqid, test_seqid)}
    reference = parse_annotation(reference_path, scored_lengths, protein_coding_only=True)
    complete_all = [tx for tx in reference["transcripts"].values()
                    if tx["CDS"] and not tx["partial"] and not reference["genes"][tx["gene_id"]]["partial"]]
    views = {"primary_transcript": primary_transcripts(reference), "all_isoform_union": complete_all}
    result = {
        "cache_path": str(cache_path.relative_to(root)), "feature_order": feature_names,
        "splits": splits, "validation_seqid": val_seqid, "test_seqid": test_seqid,
        "tile_bp": 6000, "views": {},
    }
    val_features = data[f"feat::{val_seqid}"]
    test_features = data[f"feat::{test_seqid}"]
    if val_features.shape != (lengths[val_seqid], len(feature_names)):
        raise ValueError(f"{cache_path}: validation feature shape {val_features.shape}")
    if test_features.shape != (lengths[test_seqid], len(feature_names)):
        raise ValueError(f"{cache_path}: test feature shape {test_features.shape}")
    for view_name, transcripts in views.items():
        by_seqid = defaultdict(list)
        for tx in transcripts:
            by_seqid[tx["seqid"]].append(tx)
        result["views"][view_name] = {}
        for element in ("exon", "intron", "splice_donor", "splice_acceptor"):
            feature_index = feature_names.index(element)
            val_scores = val_features[:, feature_index]
            test_scores = test_features[:, feature_index]
            val_truth = truth_mask(lengths[val_seqid], by_seqid[val_seqid], element)
            test_truth = truth_mask(lengths[test_seqid], by_seqid[test_seqid], element)
            _val_ap, threshold, val_best_f1 = score_histogram(val_scores, val_truth)
            test_ap, _unused, _unused_f1 = score_histogram(test_scores, test_truth)
            result["views"][view_name][element] = {
                "test_AUCPR": test_ap, "validation_selected_threshold": threshold,
                "validation_best_F1": val_best_f1,
                "test_F1_at_validation_threshold": score_at_threshold(test_scores, test_truth, threshold),
                "test_prevalence": float(test_truth.mean()),
            }
    return result


def render_report(metrics):
    lines = ["# M24 direct structure diagnostic", "", "## Scope audit", ""]
    for species, audit in metrics["scope_audit"].items():
        lines.append(f"- `{species}`: {audit['genome_bases']:,} bp; {audit['reference_primary_transcripts']:,} complete primary transcripts; {audit['reference_CDS_intervals']:,} CDS intervals; {audit['excluded_partial_genes']:,} partial genes excluded from exact metrics.")
    lines.extend(["", "## Structural metrics", "",
                  "Candidate rows are coordinate-only diagnostics. Their strand/phase/transcript/gene fields are `not_applicable`.", "",
                  "| method | species | exact CDS F1 | pseudo/exact CDS-chain F1 | gbF1 | intergenic FPR | gene-count ratio |",
                  "|---|---|---:|---:|---:|---:|---:|"])
    for method, species_rows in metrics["methods"].items():
        for species, row in species_rows.items():
            chain = row["structure"]["pseudo_CDS_chain"]["f1"]
            if isinstance(row["structure"].get("exact_CDS_chain"), dict):
                chain = row["structure"]["exact_CDS_chain"]["f1"]
            gb = row["gene_body"]
            lines.append(f"| {method} | {species} | {row['structure']['exact_CDS_run']['f1']:.4f} | {chain:.4f} | {gb['gene_body_F1']:.4f} | {gb['intergenic_FPR']:.4f} | {gb['predicted_gene_count_ratio']:.3f} |")
    lines.extend(["", "## SegmentNT released feature cache", "",
                  "Thresholds were selected independently on each species' validation seqid and applied once to its test seqid. The cache uses independent 6,000-bp tiles.", "",
                  "| species | view | element | test AUCPR | test F1 | prevalence |",
                  "|---|---|---|---:|---:|---:|"])
    for species, segment in metrics["segmentnt"].items():
        for view, elements in segment["views"].items():
            for element, row in elements.items():
                lines.append(f"| {species} | {view} | {element} | {row['test_AUCPR']:.4f} | {row['test_F1_at_validation_threshold']:.4f} | {row['test_prevalence']:.6f} |")
    lines.extend(["", "## Interpretation boundary", "",
                  "This report measures saved artifacts only. Candidate `+` strand and phase `0` placeholders are not ranked as structural predictions. A weak SegmentNT row applies only to the existing 6-kb tiled cache, not to a longer-context extraction or the checkpoint in general.", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--out-dir", default="outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC")
    args = parser.parse_args()
    root = Path(args.repo_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = {"experiment_id": "M24-DIRECT-STRUCTURE-DIAGNOSTIC", "scope_audit": {}, "methods": {}, "segmentnt": {}}
    references = {}
    for species in SPECIES:
        subset = root / "outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/eval_subsets" / species
        fasta = subset / "genome.fa"
        reference_path = subset / "reference.gff3"
        lengths = gene_body.fasta_lengths(fasta)
        reference = parse_annotation(reference_path, lengths, protein_coding_only=True)
        primary = primary_transcripts(reference)
        references[species] = (reference, lengths)
        metrics["scope_audit"][species] = {
            "seqids": lengths, "genome_bases": sum(lengths.values()),
            "reference_primary_transcripts": len(primary),
            "reference_CDS_intervals": len(intervals(primary, strand=False)),
            "excluded_partial_genes": sum(gene["partial"] for gene in reference["genes"].values()),
        }

    for method, (relative_dir, suffix, capability) in METHODS.items():
        metrics["methods"][method] = {}
        for species in SPECIES:
            reference, lengths = references[species]
            prediction_path = root / relative_dir / f"{species}{suffix}"
            prediction = parse_annotation(prediction_path, lengths)
            row = {
                "prediction_path": str(prediction_path.relative_to(root)),
                "retained_prediction_records": prediction["retained_records"],
                "out_of_bounds": prediction["out_of_bounds"],
                "scope": metrics["scope_audit"][species],
                "structure": exact_structure_metrics(reference, prediction, capability),
                "gene_body": gene_body_metrics(reference, prediction, lengths),
            }
            metrics["methods"][method][species] = row

    for species in SPECIES:
        metrics["segmentnt"][species] = segmentnt_metrics(root, species)

    metrics_path = out_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    report = render_report(metrics)
    (out_dir / "report.md").write_text(report)
    tracked_report = root / "reports/M24-DIRECT-STRUCTURE-DIAGNOSTIC/report.md"
    tracked_report.parent.mkdir(parents=True, exist_ok=True)
    tracked_report.write_text(report)


if __name__ == "__main__":
    main()
