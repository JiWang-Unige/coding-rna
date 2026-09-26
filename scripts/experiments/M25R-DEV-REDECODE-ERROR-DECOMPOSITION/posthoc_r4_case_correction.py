#!/usr/bin/env python3
"""One CPU-only correction of R4 diagnostic case semantics; never run inference."""

import csv
import json
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import redecode_error_decomposition as diag

ROOT = Path(__file__).resolve().parents[3]
SOURCE_ID = "M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4"
OUTPUT_ID = SOURCE_ID + "-CASE-CORRECTION"
NOT_RECOVERABLE = "NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS"
VALIDATION = {("arabidopsis_thaliana", "NC_003074.8"), ("oryza_sativa", "NC_089041.1")}


def read_json(path):
    return json.loads(path.read_text())


def source_stats(source):
    # Direct metadata comparison detects accidental source writes; no hash layer.
    return {str(p.relative_to(source)): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in source.rglob("*") if p.is_file()}


def reconstruct_states(length, traces):
    states = np.full(length, diag.m25.I, dtype=np.int8)
    previous_end = -1
    for trace in sorted(traces, key=lambda t: t["block_span"]):
        start, end = trace["block_span"]
        if not 0 <= start < end <= length or start <= previous_end:
            raise AssertionError("overlapping, adjacent or out-of-range lineage blocks")
        states[start:end] = diag.m25.G
        previous_run_end = start
        for left, right in trace["runs"]:
            if not start <= left < right <= end or left < previous_run_end:
                raise AssertionError("invalid ordered CDS runs")
            states[left:right] = diag.m25.C
            previous_run_end = right
        previous_end = end
    return states


def reachability_index(states, sequence):
    """Same transition/event/motif semantics, indexed once per chromosome strand."""
    anchors = {name: [] for name in diag.m25.BOUNDARY_NAMES}
    motifs = {name: set() for name in diag.m25.BOUNDARY_NAMES}
    for position in np.flatnonzero(states[1:] != states[:-1]) + 1:
        pair = (int(states[position - 1]), int(states[position]))
        for name in diag.m25.TRANSITIONS.get(pair, ()):
            anchor = diag.event_anchor(int(position), name)
            anchors[name].append(anchor)
            motifs[name].update(diag.motif_positions(sequence, name, anchor))
    return anchors, motifs


def indexed_reachability(reference, index):
    anchors, motifs = index
    transition = motif = True
    for name, positions in reference["events"].items():
        candidates = anchors[name]
        # Restrict only to anchors that could match a truth event, preserving order.
        nearby = sorted({a for p in positions
                         for a in candidates[bisect_left(candidates, p - diag.RADIUS_BP):
                                             bisect_right(candidates, p + diag.RADIUS_BP)]})
        matched = diag.event_match_count(positions, nearby, diag.RADIUS_BP) == len(positions)
        transition &= matched
        motif &= matched and all(p in motifs[name] for p in positions)
    return bool(transition), bool(motif)


def correct_epoch(epoch, source, output, species, references, lengths, records):
    src = source / f"epoch_{epoch}"
    old = read_json(src / "diagnostic.json")
    dst = output / f"epoch_{epoch}"
    dst.mkdir()
    sequences = {sid: species[sp]["seqs"][sid] for sp, sid in lengths}
    flat_lengths = {sid: length for (_sp, sid), length in lengths.items()}
    gff = src / "replayed_predictions.gff3"
    expected = old["candidate_lineage_reconciliation"]["final_GFF3_models"]
    validity = diag.audit_gff3(gff, sequences, flat_lengths, expected)
    if validity["valid_transcripts"] != expected or validity["audit_coverage"] != 1:
        raise AssertionError("corrected validity differs from verified all-valid result")
    with (dst / "structural_validity_corrected.jsonl").open("w") as handle:
        for row in validity.pop("transcript_ledger"):
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    parsed = diag.primary_transcripts(diag.parse_annotation(str(gff), flat_lengths, protein_coding_only=False))
    predictions = {key: [] for key in references}
    for tx in parsed:
        key = next(key for key in references if key[1] == tx["seqid"])
        predictions[key].append({"strand": tx["strand"], "cds": [(a, b) for a, b, _p in tx["CDS"]]})
    metrics = diag.m25._validation_metrics(predictions, references, lengths)
    coordinate_metrics = {name: metrics[name] for name in old["reproduction"]}
    if any(value != old["reproduction"][name]["replayed"] for name, value in coordinate_metrics.items()):
        raise AssertionError("coordinate metric changed when parsing existing GFF3")
    pred_chains = {(sp, sid, model["strand"], tuple(model["cds"]))
                   for (sp, sid), models in predictions.items() for model in models}

    with (src / "reference_attrition.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row_by_key = {(r["species"], r["seqid"], r["strand"], r["transcript_id"]): r for r in rows}
    if len(rows) != 6450 or len(row_by_key) != 6450 or set(row_by_key) != {r["key"] for r in records}:
        raise AssertionError("reference universe does not reconcile")
    assignments = {key: row["earliest_stage"] for key, row in row_by_key.items()}
    stage_counts = {s: Counter(assignments.values())[s] for s in diag.STAGES}
    if stage_counts != old["reference_attrition"]["stage_counts"]:
        raise AssertionError("source stage counts do not reconcile")
    traces = defaultdict(list)
    with (src / "candidate_lineages.jsonl").open() as handle:
        for line in handle:
            trace = json.loads(line)
            key = (trace["species"], trace["seqid"], trace["strand"])
            if key[:2] not in VALIDATION or key[2] not in {"+", "-"}:
                raise AssertionError("lineage outside approved validation scope")
            traces[key].append(trace)
    if sum(map(len, traces.values())) != old["candidate_lineage_reconciliation"]["non_intergenic_blocks"]:
        raise AssertionError("lineage count mismatch")
    grouped = defaultdict(list)
    for record in records:
        grouped[record["key"][:3]].append(record)
    for key, group in grouped.items():
        seq = sequences[key[1]]
        if key[2] == "-":
            seq = diag.m25.reverse_complement(seq)
        states = reconstruct_states(len(seq), traces[key])
        index = reachability_index(states, seq)
        for ref in group:
            transition, motif = indexed_reachability(ref, index)
            row = row_by_key[ref["key"]]
            if transition != bool(int(row["transition_reachable"])):
                raise AssertionError(f"transition mismatch: {ref['key']}")
            cds = ref["cds"]
            if key[2] == "-":
                cds = sorted((len(seq) - b, len(seq) - a) for a, b in cds)
            exact = (*key, tuple(cds)) in pred_chains
            if exact != (row["earliest_stage"] == "emitted_exact_chain"):
                raise AssertionError(f"exact identity does not reconcile: {ref['key']}")
            if exact and not (transition and motif and ref["canonical"]):
                raise AssertionError(f"exact chain outside corrected reachability: {ref['key']}")
            ref["transition_reachable"] = transition
            ref["motif_reachable"] = motif
        del states, index
        print(f"epoch={epoch} corrected {key}", flush=True)

    strata = diag.count_by_strata(records, assignments)
    for label, counts in strata.items():
        if label not in {"canonical", "noncanonical"} and counts != old["reference_attrition"]["strata"][label]:
            raise AssertionError("non-canonical stratum changed")
    canonical = [r for r in records if r["canonical"]]
    reachability = {}
    for name in ("transition_reachable", "motif_reachable"):
        count = sum(r[name] for r in records)
        reachability[name] = {"chains": count, "fraction_of_R_all": count / len(records),
                              "canonical_chains": sum(r[name] for r in canonical),
                              "fraction_of_R_canonical": sum(r[name] for r in canonical) / len(canonical)}
    if reachability["transition_reachable"]["chains"] != old["post_hoc_upper_bounds"]["transition_reachable"]["chains"]:
        raise AssertionError("transition total mismatch")
    if not stage_counts["emitted_exact_chain"] <= reachability["motif_reachable"]["chains"] <= reachability["transition_reachable"]["chains"]:
        raise AssertionError("conditional recall ceiling ordering failed")
    fields = list(rows[0])
    with (dst / "reference_case_correction.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for ref in records:
            row = dict(row_by_key[ref["key"]])
            row.update(canonical=int(ref["canonical"]), motif_reachable=int(ref["motif_reachable"]),
                       truth_assisted_exact_chain=NOT_RECOVERABLE)
            writer.writerow(row)
    result = {"epoch": epoch, "coordinate_metrics": coordinate_metrics, "coordinate_metrics_unchanged": True,
              "main_stage_counts": stage_counts, "main_stage_counts_unchanged": True,
              "validity": validity, "canonical_reference_count": len(canonical),
              "unsupported_by_frozen_canonical_decoder": len(records) - len(canonical),
              "strata": strata, "conditional_reachability": reachability,
              "truth_assisted_exact_chain": NOT_RECOVERABLE,
              "reference_total": len(records), "all_exact_references_covered": True}
    diag.save_json_atomic(dst / "correction_summary.json", result)
    return result


def main():
    source, output = ROOT / "outputs" / SOURCE_ID, ROOT / "outputs" / OUTPUT_ID
    # Slurm prepares logs/JOBID/STATUS; never reuse a previous correction payload.
    if (output / "posthoc_case_correction_summary.json").exists() or any(output.glob("epoch_*")):
        raise FileExistsError("correction output already contains result payloads")
    output.mkdir(exist_ok=True)
    before = source_stats(source)
    original = read_json(source / "stage1_diagnostic.json")
    resolved = read_json(source / "resolved_inputs.json")
    config_path = ROOT / "configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml"
    config = yaml.safe_load(config_path.read_text())
    species = diag.load_species(ROOT, config)
    references, lengths = diag.validation_truth(species)
    if set(references) != VALIDATION:
        raise AssertionError("validation chromosome scope changed")
    records = diag.build_reference_records(species, references)
    if len(records) != resolved["reference_chains"] or len(records) != 6450:
        raise AssertionError("reference count changed")
    epochs = [correct_epoch(e, source, output, species, references, lengths, records) for e in (1, 2, 3)]
    if before != source_stats(source):
        raise AssertionError("original R4 files changed during correction")
    summary = {"status": "POSTHOC_CORRECTION_COMPLETED_REVIEW_REQUIRED", "source": str(source),
               "original_invalid_structure_conclusion": "RETRACTED_CASE_FALSE_NEGATIVE",
               "original_scientific_status": original["scientific_status"],
               "m25r_frozen_contract_verdict": "DEVELOPMENT_NO_GO_UNCHANGED",
               "truth_assisted_exact_chain": NOT_RECOVERABLE, "epochs": epochs,
               "original_R4_metadata_unchanged": True, "setaria_files_read": False,
               "checkpoints_loaded": False, "model_forward_called": False,
               "threshold_or_decoder_search": False, "next_action": "STOP_FOR_REVIEW"}
    diag.save_json_atomic(output / "posthoc_case_correction_summary.json", summary)
    (output / "STATUS").write_text(summary["status"] + "\n")
    print(summary["status"], flush=True)


if __name__ == "__main__":
    main()
