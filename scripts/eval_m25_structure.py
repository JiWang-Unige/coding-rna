#!/usr/bin/env python3
"""Evaluate the frozen M25 full and unchanged-input ablation predictions."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import eval_structure_diagnostic as m24


MARKER_KEYS = {
    "config_resolved": None,
    "primary_chromosome_allowlist": "allowlist",
    "checkpoint": None,
    "validation_decode_selection": None,
    "genome_fasta": "genome_fasta",
    "full_prediction_gff3": "full_gff3",
    "ablation_prediction_gff3": "ablation_gff3",
}
WINDOW_BP = 6144
SNAP_BP = 6


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_embargo(marker_path, cli_paths):
    marker = json.loads(marker_path.read_text())
    if set(marker) != set(MARKER_KEYS):
        raise ValueError(f"{marker_path}: expected exactly {sorted(MARKER_KEYS)}")
    verified = {}
    for key, cli_name in MARKER_KEYS.items():
        record = marker[key]
        if set(record) != {"path", "sha256"}:
            raise ValueError(f"{marker_path}: {key} must contain only path and sha256")
        path = Path(record["path"]).expanduser().resolve()
        if cli_name and path != cli_paths[cli_name]:
            raise ValueError(f"{marker_path}: {key} path does not match CLI input")
        actual = sha256(path)
        if actual != record["sha256"]:
            raise ValueError(f"{path}: sha256 mismatch for {key}")
        verified[key] = {"path": str(path), "sha256": actual}
    return verified


def prediction_placeholder(annotation):
    transcripts = m24.primary_transcripts(annotation)
    bad_strand = sum(tx["strand"] not in {"+", "-"} for tx in transcripts)
    bad_phase = sum(phase not in {"0", "1", "2"} for tx in transcripts for _start, _end, phase in tx["CDS"])
    return {"bad_strand_transcripts": bad_strand, "bad_phase_CDS_intervals": bad_phase,
            "placeholder_remains": bool(bad_strand or bad_phase)}


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path):
    sequences = {}
    seqid = None
    chunks = []
    for line in path.open():
        line = line.strip()
        if line.startswith(">"):
            if seqid is not None:
                sequences[seqid] = "".join(chunks).upper()
            seqid = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if seqid is not None:
        sequences[seqid] = "".join(chunks).upper()
    return sequences


def codon_features(path):
    result = defaultdict(lambda: {"start_codon": [], "stop_codon": []})
    for line_no, line in enumerate(path.open(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"{path}:{line_no}: expected 9 columns")
        feature = fields[2]
        if feature not in {"start_codon", "stop_codon"}:
            continue
        attrs = m24.parse_attrs(fields[8])
        parents = attrs.get("Parent", "").split(",") if attrs.get("Parent") else []
        if len(parents) != 1:
            raise ValueError(f"{path}:{line_no}: codon must have one Parent")
        result[parents[0]][feature].append((fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
    return result


def ordered_cds(tx):
    return tx["CDS"] if tx["strand"] == "+" else list(reversed(tx["CDS"]))


def transcript_sequence(tx, sequences):
    parts = [sequences[tx["seqid"]][start:end] for start, end, _phase in ordered_cds(tx)]
    if tx["strand"] == "+":
        return "".join(parts)
    return "".join(reverse_complement(part) for part in parts)


def canonical_introns(tx, sequences):
    motifs = []
    for start, end in m24.introns(tx):
        intron = sequences[tx["seqid"]][start:end]
        if tx["strand"] == "-":
            intron = reverse_complement(intron)
        motifs.append((intron[:2], intron[-2:]))
    return motifs


def expected_codon_records(tx):
    first = ordered_cds(tx)[0]
    last = ordered_cds(tx)[-1]
    if tx["strand"] == "+":
        start = (tx["seqid"], first[0], first[0] + 3, "+")
        stop = (tx["seqid"], last[1] - 3, last[1], "+")
    else:
        start = (tx["seqid"], first[1] - 3, first[1], "-")
        stop = (tx["seqid"], last[0], last[0] + 3, "-")
    return start, stop


def structurally_valid(tx, sequences, codons):
    parts = ordered_cds(tx)
    if any(phase not in {"0", "1", "2"} for _start, _end, phase in parts):
        return False
    phases = [int(phase) for _start, _end, phase in parts]
    cumulative = 0
    phase_ok = True
    for (start, end, _phase), phase in zip(parts, phases):
        phase_ok &= phase == ((3 - cumulative % 3) % 3)
        cumulative += end - start
    coding = transcript_sequence(tx, sequences)
    start_record, stop_record = expected_codon_records(tx)
    codon_row = codons.get(tx["id"], {})
    codon_ok = codon_row.get("start_codon") == [start_record] and codon_row.get("stop_codon") == [stop_record]
    motif_ok = all(donor == "GT" and acceptor == "AG" for donor, acceptor in canonical_introns(tx, sequences))
    sequence_ok = len(coding) >= 6 and len(coding) % 3 == 0 and coding[:3] == "ATG" and coding[-3:] in {"TAA", "TAG", "TGA"}
    internal_stop = any(coding[index:index + 3] in {"TAA", "TAG", "TGA"}
                        for index in range(3, len(coding) - 3, 3))
    return phase_ok and codon_ok and motif_ok and sequence_ok and not internal_stop


def validity_metrics(annotation, sequences, codons):
    transcripts = m24.primary_transcripts(annotation)
    valid = sum(structurally_valid(tx, sequences, codons) for tx in transcripts)
    return {"valid_complete_transcripts": valid, "complete_transcripts": len(transcripts),
            "fraction": valid / len(transcripts) if transcripts else 0.0}


def cds_overlap(a, b):
    overlap = 0
    for a_start, a_end, _phase in a["CDS"]:
        for b_start, b_end, _phase in b["CDS"]:
            overlap += max(0, min(a_end, b_end) - max(a_start, b_start))
    return overlap


def matched_gene_strand_accuracy(reference_tx, predicted_tx):
    candidates = []
    for ref in reference_tx:
        for pred in predicted_tx:
            if ref["seqid"] == pred["seqid"]:
                overlap = cds_overlap(ref, pred)
                if overlap:
                    candidates.append((-overlap, ref["id"], pred["id"], ref, pred))
    used_ref = set()
    used_pred = set()
    matches = []
    for _negative_overlap, ref_id, pred_id, ref, pred in sorted(candidates):
        if ref_id not in used_ref and pred_id not in used_pred:
            used_ref.add(ref_id)
            used_pred.add(pred_id)
            matches.append(ref["strand"] == pred["strand"])
    return {"matched_genes": len(matches), "accuracy": sum(matches) / len(matches) if matches else m24.NA}


def phase_accuracy(reference_tx, predicted_tx):
    ref = {(tx["seqid"], start, end, tx["strand"]): phase
           for tx in reference_tx for start, end, phase in tx["CDS"]}
    pred = {(tx["seqid"], start, end, tx["strand"]): phase
            for tx in predicted_tx for start, end, phase in tx["CDS"]}
    keys = [key for key in ref.keys() & pred.keys() if ref[key] in {"0", "1", "2"} and pred[key] in {"0", "1", "2"}]
    return {"exact_matched_CDS": len(keys),
            "accuracy": sum(ref[key] == pred[key] for key in keys) / len(keys) if keys else m24.NA}


def boundary_events(transcripts):
    events = defaultdict(list)
    for tx in transcripts:
        parts = ordered_cds(tx)
        if tx["strand"] == "+":
            events[(tx["seqid"], tx["strand"], "start")].append(parts[0][0])
            events[(tx["seqid"], tx["strand"], "stop")].append(parts[-1][1])
            for left, right in zip(parts, parts[1:]):
                events[(tx["seqid"], tx["strand"], "donor")].append(left[1])
                events[(tx["seqid"], tx["strand"], "acceptor")].append(right[0])
        else:
            events[(tx["seqid"], tx["strand"], "start")].append(parts[0][1])
            events[(tx["seqid"], tx["strand"], "stop")].append(parts[-1][0])
            for upstream, downstream in zip(parts, parts[1:]):
                events[(tx["seqid"], tx["strand"], "donor")].append(upstream[0])
                events[(tx["seqid"], tx["strand"], "acceptor")].append(downstream[1])
    return events


def boundary_diagnostics(reference_tx, predicted_tx):
    reference = boundary_events(reference_tx)
    predicted = boundary_events(predicted_tx)
    result = {}
    for kind in ("start", "stop", "donor", "acceptor"):
        offsets = []
        reference_count = 0
        for key, ref_positions in reference.items():
            if key[2] != kind:
                continue
            reference_count += len(ref_positions)
            choices = predicted.get(key, [])
            for position in ref_positions:
                if choices:
                    nearest = min(choices, key=lambda value: (abs(value - position), value))
                    if abs(nearest - position) <= SNAP_BP:
                        offsets.append(nearest - position)
        counts = Counter(offsets)
        result[kind] = {
            "reference_boundaries": reference_count,
            "matched_within_6bp": len(offsets),
            "signed_offset_counts": {str(offset): counts[offset] for offset in range(-SNAP_BP, SNAP_BP + 1)},
            "modulo_6_counts": {str(residue): sum(count for offset, count in counts.items() if offset % 6 == residue)
                                for residue in range(6)},
        }
    return result


def stratum_metrics(reference_tx, predicted_tx):
    ref_intervals = m24.intervals(reference_tx)
    pred_intervals = m24.intervals(predicted_tx)
    ref_chains = m24.chains(reference_tx)
    pred_chains = m24.chains(predicted_tx)
    return {
        "reference_transcripts": len(reference_tx),
        "exact_CDS_interval_recall": len(ref_intervals & pred_intervals) / len(ref_intervals) if ref_intervals else m24.NA,
        "exact_CDS_chain_recall": len(ref_chains & pred_chains) / len(ref_chains) if ref_chains else m24.NA,
    }


def fixed_strata(reference_tx, predicted_tx, sequences):
    strata = {
        "noncanonical_splice": [tx for tx in reference_tx if any(motif != ("GT", "AG") for motif in canonical_introns(tx, sequences))],
        "span_gt_6144bp": [tx for tx in reference_tx if tx["CDS"][-1][1] - tx["CDS"][0][0] > WINDOW_BP],
        "CDS_boundary_within_6bp_of_6144bp_tile_edge": [
            tx for tx in reference_tx
            if any(min(position % WINDOW_BP, WINDOW_BP - position % WINDOW_BP) <= SNAP_BP
                   for start, end, _phase in tx["CDS"] for position in (start, end))
        ],
    }
    return {name: stratum_metrics(transcripts, predicted_tx) for name, transcripts in strata.items()}


def evaluate(reference, prediction, lengths, sequences, prediction_path):
    ref_tx = m24.primary_transcripts(reference)
    pred_tx = m24.primary_transcripts(prediction)
    ref_cds = {(seqid, start, end, strand) for seqid, start, end, strand, _phase in m24.intervals(ref_tx)}
    pred_cds = {(seqid, start, end, strand) for seqid, start, end, strand, _phase in m24.intervals(pred_tx)}
    chain = m24.prf(m24.chains(pred_tx), m24.chains(ref_tx))
    body = m24.gene_body_metrics(reference, prediction, lengths)
    return {
        "strand_aware_exact_CDS_interval": m24.prf(pred_cds, ref_cds),
        "exact_CDS_chain": chain,
        "exact_coding_transcript": chain,
        "exact_coding_gene": chain,
        "matched_gene_strand": matched_gene_strand_accuracy(ref_tx, pred_tx),
        "exact_matched_CDS_phase": phase_accuracy(ref_tx, pred_tx),
        "structural_validity": validity_metrics(prediction, sequences, codon_features(prediction_path)),
        "intergenic_FPR": body["intergenic_FPR"],
        "predicted_gene_count_ratio": body["predicted_gene_count_ratio"],
        "placeholder_audit": prediction_placeholder(prediction),
        "signed_boundary_offset_and_modulo_6": boundary_diagnostics(ref_tx, pred_tx),
        "fixed_strata": fixed_strata(ref_tx, pred_tx, sequences),
    }


def numeric_at(row, *keys):
    value = row
    for key in keys:
        value = value[key]
    return value if isinstance(value, (int, float)) else float("-inf")


def decide(full, ablation):
    gains = {
        "exact_CDS_interval_F1": numeric_at(full, "strand_aware_exact_CDS_interval", "f1") - numeric_at(ablation, "strand_aware_exact_CDS_interval", "f1"),
        "exact_CDS_chain_F1": numeric_at(full, "exact_CDS_chain", "f1") - numeric_at(ablation, "exact_CDS_chain", "f1"),
        "intergenic_FPR": numeric_at(full, "intergenic_FPR") - numeric_at(ablation, "intergenic_FPR"),
    }
    success = {
        "exact_CDS_interval_F1_gte_0.80": numeric_at(full, "strand_aware_exact_CDS_interval", "f1") >= 0.80,
        "exact_CDS_chain_F1_gte_0.55": numeric_at(full, "exact_CDS_chain", "f1") >= 0.55,
        "exact_coding_gene_F1_gte_0.50": numeric_at(full, "exact_coding_gene", "f1") >= 0.50,
        "matched_gene_strand_accuracy_gte_0.98": numeric_at(full, "matched_gene_strand", "accuracy") >= 0.98,
        "phase_accuracy_gte_0.90": numeric_at(full, "exact_matched_CDS_phase", "accuracy") >= 0.90,
        "structural_validity_gte_0.99": numeric_at(full, "structural_validity", "fraction") >= 0.99,
        "intergenic_FPR_lte_0.020": numeric_at(full, "intergenic_FPR") <= 0.020,
        "gene_count_ratio_0.80_to_1.20": 0.80 <= numeric_at(full, "predicted_gene_count_ratio") <= 1.20,
        "interval_ablation_gain_gte_0.10": gains["exact_CDS_interval_F1"] >= 0.10,
        "chain_ablation_gain_gte_0.10": gains["exact_CDS_chain_F1"] >= 0.10,
        "FPR_ablation_delta_lte_0.005": gains["intergenic_FPR"] <= 0.005,
    }
    immediate_stop = {
        "placeholder_strand_or_phase": full["placeholder_audit"]["placeholder_remains"],
        "structural_validity_lt_0.95": numeric_at(full, "structural_validity", "fraction") < 0.95,
        "exact_CDS_interval_F1_lt_0.60": numeric_at(full, "strand_aware_exact_CDS_interval", "f1") < 0.60,
        "exact_CDS_chain_F1_lt_0.30": numeric_at(full, "exact_CDS_chain", "f1") < 0.30,
        "exact_coding_gene_F1_lt_0.25": numeric_at(full, "exact_coding_gene", "f1") < 0.25,
        "intergenic_FPR_gt_0.030": numeric_at(full, "intergenic_FPR") > 0.030,
        "gene_count_ratio_outside_0.50_to_1.50": not 0.50 <= numeric_at(full, "predicted_gene_count_ratio") <= 1.50,
        "interval_ablation_gain_lt_0.10": gains["exact_CDS_interval_F1"] < 0.10,
        "chain_ablation_gain_lt_0.10": gains["exact_CDS_chain_F1"] < 0.10,
    }
    return {"full_minus_ablation": gains, "success_checks": success,
            "immediate_stop_checks": immediate_stop,
            "decision": "PASSED_DISCOVERY_GATE" if all(success.values()) and not any(immediate_stop.values()) else "STOP_M25_BRANCH"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--reference-gff3", required=True)
    parser.add_argument("--full-gff3", required=True)
    parser.add_argument("--ablation-gff3", required=True)
    parser.add_argument("--embargo-marker", required=True)
    parser.add_argument("--allowlist", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()
    paths = {name: Path(value).expanduser().resolve() for name, value in vars(args).items()}

    verified = verify_embargo(paths["embargo_marker"], paths)
    all_sequences = read_fasta(paths["genome_fasta"])
    allowlists = json.loads(paths["allowlist"].read_text())
    setaria_seqids = allowlists["setaria_viridis"]
    missing = [seqid for seqid in setaria_seqids if seqid not in all_sequences]
    if missing:
        raise ValueError(f"primary chromosome seqids absent from Setaria FASTA: {missing}")
    sequences = {seqid: all_sequences[seqid] for seqid in setaria_seqids}
    lengths = {seqid: len(sequence) for seqid, sequence in sequences.items()}
    full_annotation = m24.parse_annotation(paths["full_gff3"], lengths)
    ablation_annotation = m24.parse_annotation(paths["ablation_gff3"], lengths)
    reference = m24.parse_annotation(paths["reference_gff3"], lengths, protein_coding_only=True)

    full = evaluate(reference, full_annotation, lengths, sequences, paths["full_gff3"])
    ablation = evaluate(reference, ablation_annotation, lengths, sequences, paths["ablation_gff3"])
    result = {
        "experiment_id": "M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0",
        "embargo_hashes_verified_before_reference_parse": verified,
        "reference_scope": {
            "complete_primary_protein_coding_transcripts": len(m24.primary_transcripts(reference)),
            "genome_bases": sum(lengths.values()),
        },
        "full": full,
        "ablation": ablation,
        "gate": decide(full, ablation),
    }
    paths["output_json"].parent.mkdir(parents=True, exist_ok=True)
    paths["output_json"].write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
