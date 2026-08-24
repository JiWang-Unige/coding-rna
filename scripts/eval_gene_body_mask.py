#!/usr/bin/env python3
"""Evaluate gene-body masks from GTF spans with profile-aware FPR thresholds."""

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


# Gene-body span construction is profile of which features define the per-transcript
# min..max extent. Cross-tool comparability (M1, Helixer-vs-Tiberius finding):
#   transcript = CDS+exon+intron+start/stop -> full mRNA span INCL UTR (via exon). Fair only
#                between tools that all emit UTR; Tiberius (CDS-only caller) has none, so this
#                penalizes CDS-only callers and is NOT apples-to-apples across tools.
#   cds        = CDS(+start/stop) only -> coding-region span, the one layer EVERY gene caller
#                emits (RefSeq, Tiberius, Helixer, ANNEVO). Canonical fair common denominator
#                for protein-coding annotation; removes the tool-dependent UTR confound.
SPAN_FEATURES_BY_MODE = {
    "transcript": {"CDS", "exon", "intron", "start_codon", "stop_codon"},
    "cds": {"CDS", "start_codon", "stop_codon"},
}
TRANSCRIPT_FEATURES = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "rRNA", "tRNA"}


def parse_attrs(attr_text):
    attrs = {}
    for part in attr_text.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            attrs[key] = val.strip().strip('"')
        elif "=" in part:
            key, val = part.split("=", 1)
            attrs[key] = val.strip().strip('"')
    return attrs


def fasta_lengths(path):
    lengths = {}
    current = None
    total = 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = total
                current = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if current is not None:
        lengths[current] = total
    return lengths


def collect_spans(gtf_path, span_features):
    counts = Counter()
    spans = {}
    genes = set()
    transcripts = set()
    for line_no, line in enumerate(open(gtf_path), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        seqid, _source, feature, start_s, end_s, _score, strand, _phase, attrs_s = fields
        counts[feature] += 1
        attrs = parse_attrs(attrs_s)
        attr_id = attrs.get("ID")
        parent_id = attrs.get("Parent")
        gene_id = attrs.get("gene_id")
        transcript_id = attrs.get("transcript_id")
        if not gene_id and feature == "gene":
            gene_id = attr_id
        if not transcript_id and feature in TRANSCRIPT_FEATURES:
            transcript_id = attr_id
        if gene_id:
            genes.add(gene_id)
        if transcript_id:
            transcripts.add(transcript_id)
        if feature not in span_features:
            continue
        group_id = transcript_id or gene_id or parent_id or attr_id or f"line:{line_no}"
        start = int(start_s) - 1
        end = int(end_s)
        key = (seqid, strand, group_id)
        if key not in spans:
            spans[key] = [start, end]
        else:
            spans[key][0] = min(spans[key][0], start)
            spans[key][1] = max(spans[key][1], end)
    by_seqid = defaultdict(list)
    for seqid, _strand, _group_id in spans:
        start, end = spans[(seqid, _strand, _group_id)]
        by_seqid[seqid].append((start, end))
    return {
        "intervals_by_seqid": {k: merge_intervals(v) for k, v in by_seqid.items()},
        "feature_counts": dict(sorted(counts.items())),
        "gene_count": len(genes),
        "transcript_count": len(transcripts),
        "span_group_count": len(spans),
    }


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if start < 0 or end < start:
            raise ValueError(f"invalid interval: {(start, end)}")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_length(by_seqid):
    return sum(end - start for intervals in by_seqid.values() for start, end in intervals)


def intersection_length(a_by_seqid, b_by_seqid):
    total = 0
    for seqid in set(a_by_seqid) & set(b_by_seqid):
        a = a_by_seqid[seqid]
        b = b_by_seqid[seqid]
        i = j = 0
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                total += end - start
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
    return total


def f1(precision, recall):
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def load_previous(path):
    if not path:
        return {}
    with open(path) as handle:
        data = json.load(handle)
    return {k: v for k, v in data.items() if isinstance(v, (int, float, str, bool, dict))}


def finite_float(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-gtf", required=True)
    parser.add_argument("--prediction-gtf", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--profile", default="smoke", choices=["smoke", "screen", "full", "scale"])
    parser.add_argument("--span-mode", default="transcript", choices=["transcript", "cds"],
                        help="gene-body span definition: 'transcript' (CDS+exon+intron+start/stop, "
                             "incl UTR via exon) or 'cds' (CDS-only, the fair cross-tool common "
                             "layer for protein-coding annotation). Use 'cds' for cross-tool screen.")
    parser.add_argument("--previous-metrics", default=None)
    parser.add_argument("--smoke-screen-fpr-threshold", type=float, default=0.02)
    parser.add_argument("--full-scale-fpr-threshold", type=float, default=0.01)
    parser.add_argument("--sensitivity-thresholds", default="0.005,0.01,0.02")
    args = parser.parse_args()

    lengths = fasta_lengths(args.genome_fasta)
    genome_bases = sum(lengths.values())
    span_features = SPAN_FEATURES_BY_MODE[args.span_mode]
    ref = collect_spans(args.reference_gtf, span_features)
    pred = collect_spans(args.prediction_gtf, span_features)

    ref_len = interval_length(ref["intervals_by_seqid"])
    pred_len = interval_length(pred["intervals_by_seqid"])
    overlap = intersection_length(ref["intervals_by_seqid"], pred["intervals_by_seqid"])
    pred_only = pred_len - overlap

    precision = overlap / pred_len if pred_len else 0.0
    recall = overlap / ref_len if ref_len else 0.0
    gene_body_f1 = f1(precision, recall)

    # --- INTERGENIC = complement of the FULL-TRANSCRIPT span (incl UTR), INDEPENDENT of span_mode ---
    # (revise-goal 2026-06-11): a base in a UTR is GENIC, not intergenic, so the intergenic region
    # must be genome - full_transcript_span (NOT genome - CDS_span). intergenic_FP = predicted-genic
    # bases that fall OUTSIDE every reference transcript (truly intergenic). primary metric =
    # intergenic_specificity = 1 - intergenic_FPR. Old CDS-complement value kept as a diagnostic.
    ref_full = collect_spans(args.reference_gtf, SPAN_FEATURES_BY_MODE["transcript"])
    ref_full_len = interval_length(ref_full["intervals_by_seqid"])
    true_intergenic_bases = genome_bases - ref_full_len
    pred_in_transcript = intersection_length(pred["intervals_by_seqid"], ref_full["intervals_by_seqid"])
    intergenic_fp = pred_len - pred_in_transcript           # predicted-genic bases outside ALL ref transcripts
    intergenic_fpr = intergenic_fp / true_intergenic_bases if true_intergenic_bases > 0 else 0.0
    intergenic_specificity = 1.0 - intergenic_fpr
    intergenic_fpr_cds_complement_diag = (pred_only / (genome_bases - ref_len)) if (genome_bases - ref_len) > 0 else 0.0

    threshold = args.smoke_screen_fpr_threshold if args.profile in {"smoke", "screen"} else args.full_scale_fpr_threshold
    sensitivities = [float(x) for x in args.sensitivity_thresholds.split(",") if x.strip()]

    metrics = load_previous(args.previous_metrics)
    metrics.update({
        "experiment_id": args.experiment_id,
        "profile": args.profile,
        "primary_metric": "intergenic_specificity",
        "semantic_success": True,
        "span_mode": args.span_mode,
        "gene_body_mask_mode": (
            "cds_only_span_from_CDS_start_stop" if args.span_mode == "cds"
            else "symmetric_transcript_span_from_CDS_exon_intron_start_stop"
        ),
        # PRIMARY (revise-goal 2026-06-11): intergenic = genome - FULL-transcript span (incl UTR),
        # decoupled from the gene-body-F1 span_mode. specificity = 1 - intergenic_FPR (higher=better).
        "intergenic_specificity": intergenic_specificity,
        "intergenic_FPR": intergenic_fpr,
        "intergenic_definition": "complement_of_full_transcript_span_incl_UTR",
        # secondary: gene-body F1 (CDS span) — demoted from primary; still the anti-degeneracy signal.
        "intergenic_FPR_threshold_used": threshold,
        "intergenic_guardrail_pass": intergenic_fpr <= threshold,
        "constrained_gene_body_F1": gene_body_f1 if intergenic_fpr <= threshold else 0.0,
        "gene_body_F1_unconstrained": gene_body_f1,
        "gene_body_precision": precision,
        "gene_body_recall": recall,
        "genome_bases": genome_bases,
        "reference_gene_body_bases": ref_len,
        "reference_full_transcript_bases": ref_full_len,
        "predicted_gene_body_bases": pred_len,
        "gene_body_overlap_bases": overlap,
        "predicted_intergenic_false_positive_bases": intergenic_fp,
        "reference_intergenic_bases": true_intergenic_bases,
        "intergenic_FPR_cds_complement_diag": intergenic_fpr_cds_complement_diag,
        "predicted_genic_bases_in_utr_or_cds_diag": pred_only,
        "reference_feature_counts": ref["feature_counts"],
        "prediction_feature_counts": pred["feature_counts"],
        "reference_gene_count": ref["gene_count"],
        "reference_transcript_count": ref["transcript_count"],
        "prediction_transcript_count": pred["transcript_count"],
        "reference_span_group_count": ref["span_group_count"],
        "prediction_span_group_count": pred["span_group_count"],
        "predicted_gene_count": pred["gene_count"],
        "predicted_gene_count_ratio_vs_reference": (
            pred["gene_count"] / ref["gene_count"] if ref["gene_count"] else 0.0
        ),
        "predicted_transcript_count_ratio_vs_reference": (
            pred["transcript_count"] / ref["transcript_count"] if ref["transcript_count"] else 0.0
        ),
        "nucleotide_gene_body_F1_drop_vs_anchor": 0.0,
        "notes": (
            "Gene-body masks use the same transcript-span construction for reference and "
            "prediction; sensitivity thresholds are reported separately."
        ),
    })
    for thr in sensitivities:
        suffix = str(thr).rstrip("0").rstrip(".")
        metrics[f"intergenic_guardrail_pass_at_{suffix}"] = intergenic_fpr <= thr
        metrics[f"constrained_gene_body_F1_at_{suffix}"] = gene_body_f1 if intergenic_fpr <= thr else 0.0

    clean_metrics = {k: finite_float(v) for k, v in metrics.items()}
    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        json.dump(clean_metrics, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
