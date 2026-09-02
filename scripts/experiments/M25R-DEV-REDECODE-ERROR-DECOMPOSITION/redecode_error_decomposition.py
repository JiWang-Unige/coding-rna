#!/usr/bin/env python3
"""Read-only M25R development re-decode and decoder attrition diagnostic."""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
for _path in (REPO_ROOT, REPO_ROOT / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import eval_m25_structure as m25_eval  # noqa: E402
from eval_structure_diagnostic import parse_annotation, primary_transcripts  # noqa: E402
from src.foundation_probe import train_generanno_structural_heads as m25
from src.foundation_probe.train_generanno_lora_3class import _clean, _tokenize_window
from src.foundation_probe.train_probe_head import _ConvLSTMHead
from src.screen_anchor import data as screen_data


EXPERIMENT_ID = "M25R-DEV-REDECODE-ERROR-DECOMPOSITION"
EXPECTED_REFERENCE_CHAINS = 6450
EXPECTED_TRAIN_WINDOWS = 1536
WINDOW_BP = 6144
RADIUS_BP = 6
TOLERANCES = (0, 1, 3, 6)
VALIDITY_COMPONENTS = (
    "parent_linkage",
    "phase_values",
    "phase_continuity",
    "start_codon_feature",
    "stop_codon_feature",
    "splice_motif",
    "minimum_CDS_length",
    "frame_length",
    "start_codon_sequence",
    "stop_codon_sequence",
    "internal_stop_absent",
)
STAGES = (
    "region_state_path",
    "non_intergenic_block",
    "ordered_CDS_runs",
    "legal_terminal_transitions",
    "start_stop_motif_candidates",
    "ordered_donor_acceptor_candidates",
    "learned_boundary_choice",
    "phase_check",
    "complete_ORF_internal_stop_check",
    "boundary_threshold_filter",
    "emitted_exact_chain",
)


def save_json_atomic(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def select_epoch_rows(grid):
    rows_by_epoch = defaultdict(list)
    for row in grid["rows"]:
        rows_by_epoch[int(row["epoch"])].append(row)
    selected = {}
    for epoch in (1, 2, 3):
        rows = rows_by_epoch[epoch]
        selected[epoch] = max(rows, key=lambda row: tuple(row["rank"]))
    return selected


def motif_positions(sequence, event, anchor, radius=RADIUS_BP):
    motifs = m25.MOTIFS[event]
    width = len(next(iter(motifs)))
    return [
        position
        for position in range(max(0, anchor - radius), min(len(sequence) - width, anchor + radius) + 1)
        if sequence[position:position + width] in motifs
    ]


def non_intergenic_blocks(states):
    blocks = []
    position = 0
    while position < len(states):
        while position < len(states) and states[position] == m25.I:
            position += 1
        block_start = position
        while position < len(states) and states[position] != m25.I:
            position += 1
        if block_start == position:
            continue
        runs = []
        cursor = block_start
        while cursor < position:
            while cursor < position and states[cursor] != m25.C:
                cursor += 1
            start = cursor
            while cursor < position and states[cursor] == m25.C:
                cursor += 1
            if start < cursor:
                runs.append((start, cursor))
        blocks.append({"span": (block_start, position), "runs": runs})
    return blocks


def model_key(model):
    return (
        tuple(tuple(item) for item in model["cds"]),
        tuple(model["phase"]),
        tuple(model["start_codon"]),
        tuple(model["stop_codon"]),
        tuple((name, float(model["boundary_scores"][name])) for name in m25.BOUNDARY_NAMES),
    )


def trace_decode_orientation(sequence, region_logits, boundary_logits, phase_logits, thresholds,
                             reverse_mapped=False):
    """Mirror the frozen validation decoder and retain every block lineage."""
    sequence = sequence.upper()
    states = m25.region_state_path(region_logits, float(thresholds["region"]))
    boundary_probability = m25._sigmoid(boundary_logits)
    phase_class = np.asarray(phase_logits).argmax(axis=-1)
    traces = []
    prefilter_models = []
    emitted_models = []

    for serial, block in enumerate(non_intergenic_blocks(states), 1):
        runs = block["runs"]
        trace = {
            "lineage_id": serial,
            "block_span": block["span"],
            "runs": runs,
            "failure_stage": None,
            "candidate_positions": {name: [] for name in m25.BOUNDARY_NAMES},
            "chosen_positions": {name: [] for name in m25.BOUNDARY_NAMES},
            "phase_checks": [],
        }
        traces.append(trace)
        if not runs:
            trace["failure_stage"] = "ordered_CDS_runs"
            continue

        first_start, final_end = runs[0][0], runs[-1][1]
        if first_start == 0 or final_end == len(states):
            trace["failure_stage"] = "legal_terminal_transitions"
            continue
        if "start" not in m25.TRANSITIONS.get((int(states[first_start - 1]), int(states[first_start])), ()):
            trace["failure_stage"] = "legal_terminal_transitions"
            continue
        if "stop" not in m25.TRANSITIONS.get((int(states[final_end - 1]), int(states[final_end])), ()):
            trace["failure_stage"] = "legal_terminal_transitions"
            continue

        anchors = {
            "start": [first_start],
            "stop": [final_end - 3],
            "donor": [left[1] for left, _right in zip(runs, runs[1:])],
            "acceptor": [right[0] - 2 for _left, right in zip(runs, runs[1:])],
        }
        for name in m25.BOUNDARY_NAMES:
            trace["candidate_positions"][name] = [motif_positions(sequence, name, anchor) for anchor in anchors[name]]

        start_pick = m25._select_motif(
            sequence, "start", first_start, boundary_probability[:, 0], RADIUS_BP, False, reverse_mapped
        )
        stop_pick = m25._select_motif(
            sequence, "stop", final_end - 3, boundary_probability[:, 1], RADIUS_BP, False, reverse_mapped
        )
        if start_pick is None or stop_pick is None:
            trace["failure_stage"] = "start_stop_motif_candidates"
            continue

        chosen = {"start": [start_pick], "stop": [stop_pick], "donor": [], "acceptor": []}
        failed = False
        for left, right in zip(runs, runs[1:]):
            donor = m25._select_motif(
                sequence, "donor", left[1], boundary_probability[:, 2], RADIUS_BP, False, reverse_mapped
            )
            acceptor = m25._select_motif(
                sequence, "acceptor", right[0] - 2, boundary_probability[:, 3], RADIUS_BP, False, reverse_mapped
            )
            if donor is None or acceptor is None or donor[0] + 2 > acceptor[0]:
                trace["failure_stage"] = "ordered_donor_acceptor_candidates"
                failed = True
                break
            chosen["donor"].append(donor)
            chosen["acceptor"].append(acceptor)
        if failed:
            continue
        trace["chosen_positions"] = {
            name: [int(position) for position, _score in picks] for name, picks in chosen.items()
        }

        cds = []
        for index in range(len(runs)):
            start = start_pick[0] if index == 0 else chosen["acceptor"][index - 1][0] + 2
            end = stop_pick[0] + 3 if index == len(runs) - 1 else chosen["donor"][index][0]
            if end <= start:
                trace["failure_stage"] = "ordered_donor_acceptor_candidates"
                failed = True
                break
            cds.append((start, end))
        if failed:
            continue

        phases = []
        coding_length = 0
        for index, (start, end) in enumerate(cds):
            expected = 0 if index == 0 else (3 - coding_length % 3) % 3
            predicted = int(phase_class[start])
            trace["phase_checks"].append({
                "position": int(start), "expected": int(expected + 1), "predicted": predicted,
            })
            phases.append(expected)
            if predicted != expected + 1:
                trace["failure_stage"] = "phase_check"
                failed = True
                break
            coding_length += end - start
        if failed:
            continue

        coding_sequence = "".join(sequence[start:end] for start, end in cds)
        if (
            len(coding_sequence) < 6
            or len(coding_sequence) % 3
            or not coding_sequence.startswith("ATG")
            or coding_sequence[-3:] not in m25.STOP_CODONS
            or any(
                coding_sequence[offset:offset + 3] in m25.STOP_CODONS
                for offset in range(3, len(coding_sequence) - 3, 3)
            )
        ):
            trace["failure_stage"] = "complete_ORF_internal_stop_check"
            continue

        boundary_scores = {
            name: min((score for _position, score in picks), default=1.0)
            for name, picks in chosen.items()
        }
        model = {
            "cds": cds,
            "phase": phases,
            "start_codon": (start_pick[0], start_pick[0] + 3),
            "stop_codon": (stop_pick[0], stop_pick[0] + 3),
            "boundary_scores": boundary_scores,
            "region_span": (runs[0][0], runs[-1][1]),
        }
        trace["prefilter_model"] = model
        prefilter_models.append(model)
        if any(boundary_scores[name] < float(thresholds[name]) for name in m25.BOUNDARY_NAMES):
            trace["failure_stage"] = "boundary_threshold_filter"
            continue
        trace["emitted_model"] = model
        emitted_models.append(model)

    production_prefilter = m25.decode_orientation(
        sequence,
        region_logits,
        boundary_logits,
        phase_logits,
        {"region": thresholds["region"], **{name: 0.0 for name in m25.BOUNDARY_NAMES}},
        reverse_mapped=reverse_mapped,
    )
    production_emitted = m25._filter_models(production_prefilter, thresholds)
    if sorted(map(model_key, prefilter_models)) != sorted(map(model_key, production_prefilter)):
        raise AssertionError("trace pre-filter models do not reproduce the frozen decoder")
    if sorted(map(model_key, emitted_models)) != sorted(map(model_key, production_emitted)):
        raise AssertionError("trace emitted models do not reproduce the frozen boundary filter")
    return states, traces, prefilter_models, emitted_models


def event_match_count(reference, predicted, tolerance):
    reference = np.sort(np.asarray(reference, dtype=np.int64))
    predicted = np.sort(np.asarray(predicted, dtype=np.int64))
    j = matched = 0
    for position in reference:
        j = max(j, int(np.searchsorted(predicted, position - tolerance, side="left")))
        if j < len(predicted) and predicted[j] <= position + tolerance:
            matched += 1
            j += 1
    return matched


def float16_order_key(values):
    bits = np.asarray(values, dtype=np.float16).view(np.uint16)
    return np.where(bits & 0x8000, np.bitwise_not(bits), bits ^ 0x8000).astype(np.uint16)


def new_raw_accumulator():
    return {
        "region_confusion": np.zeros((3, 3), dtype=np.int64),
        "boundary_total_hist": np.zeros((4, 65536), dtype=np.int64),
        "boundary_positive_hist": np.zeros((4, 65536), dtype=np.int64),
        "event_reference": np.zeros(4, dtype=np.int64),
        "event_predicted": np.zeros(4, dtype=np.int64),
        "event_matched": np.zeros((4, len(TOLERANCES)), dtype=np.int64),
        "phase_correct": 0,
        "phase_total": 0,
        "CDS_start_phase_correct": 0,
        "CDS_start_phase_total": 0,
        "bases": 0,
        "examples": 0,
    }


def add_raw_metrics(accumulator, region_logits, boundary_logits, phase_logits, targets, thresholds):
    region_truth = np.asarray(targets["region"])
    boundary_truth = np.asarray(targets["boundary"])
    phase_truth = np.asarray(targets["phase"])
    structural_mask = np.asarray(targets["structural_mask"], dtype=bool)
    region_predicted = np.asarray(region_logits).argmax(axis=-1)
    accumulator["region_confusion"] += np.bincount(
        region_truth * 3 + region_predicted, minlength=9
    ).reshape(3, 3)
    accumulator["bases"] += int(region_truth.size)
    accumulator["examples"] += 1

    boundary_probability = m25._sigmoid(boundary_logits)
    for index, name in enumerate(m25.BOUNDARY_NAMES):
        scores = np.asarray(boundary_logits[:, index], dtype=np.float16)[structural_mask]
        labels = boundary_truth[:, index][structural_mask].astype(np.int64)
        keys = float16_order_key(scores).astype(np.int64)
        accumulator["boundary_total_hist"][index] += np.bincount(keys, minlength=65536)
        accumulator["boundary_positive_hist"][index] += np.bincount(
            keys, weights=labels, minlength=65536
        ).astype(np.int64)
        reference = np.flatnonzero((boundary_truth[:, index] > 0) & structural_mask)
        predicted = np.flatnonzero((boundary_probability[:, index] >= float(thresholds[name])) & structural_mask)
        accumulator["event_reference"][index] += len(reference)
        accumulator["event_predicted"][index] += len(predicted)
        for tolerance_index, tolerance in enumerate(TOLERANCES):
            accumulator["event_matched"][index, tolerance_index] += event_match_count(
                reference, predicted, tolerance
            )

    phase_predicted = np.asarray(phase_logits).argmax(axis=-1)
    phase_mask = structural_mask & (phase_truth > 0)
    accumulator["phase_correct"] += int((phase_predicted[phase_mask] == phase_truth[phase_mask]).sum())
    accumulator["phase_total"] += int(phase_mask.sum())
    cds_starts = list(np.flatnonzero((boundary_truth[:, 0] > 0) & structural_mask))
    cds_starts.extend(int(position + 2) for position in np.flatnonzero((boundary_truth[:, 3] > 0) & structural_mask)
                      if position + 2 < len(phase_truth) and structural_mask[position + 2])
    accumulator["CDS_start_phase_correct"] += sum(
        phase_predicted[position] == phase_truth[position] for position in cds_starts if phase_truth[position] > 0
    )
    accumulator["CDS_start_phase_total"] += sum(phase_truth[position] > 0 for position in cds_starts)


def average_precision_from_hist(total, positive):
    positives = int(positive.sum())
    if not positives:
        return "not_applicable"
    total = total[::-1]
    positive = positive[::-1]
    cumulative_total = np.cumsum(total)
    cumulative_positive = np.cumsum(positive)
    occupied = total > 0
    precision = cumulative_positive[occupied] / cumulative_total[occupied]
    return float(np.sum(precision * positive[occupied] / positives))


def safe_fraction(numerator, denominator):
    return float(numerator / denominator) if denominator else "not_applicable"


def finalize_raw_metrics(accumulator):
    confusion = accumulator["region_confusion"]
    classes = {}
    f1_values = []
    for index, name in enumerate(("intergenic", "CDS", "gene_body_non_CDS")):
        tp = int(confusion[index, index])
        predicted = int(confusion[:, index].sum())
        reference = int(confusion[index, :].sum())
        precision = safe_fraction(tp, predicted)
        recall = safe_fraction(tp, reference)
        f1 = (2 * precision * recall / (precision + recall)
              if isinstance(precision, float) and isinstance(recall, float) and precision + recall else 0.0)
        classes[name] = {"precision": precision, "recall": recall, "f1": f1, "reference_bases": reference}
        f1_values.append(f1)
    boundary = {}
    for index, name in enumerate(m25.BOUNDARY_NAMES):
        reference = int(accumulator["event_reference"][index])
        predicted = int(accumulator["event_predicted"][index])
        boundary[name] = {
            "AUCPR": average_precision_from_hist(
                accumulator["boundary_total_hist"][index], accumulator["boundary_positive_hist"][index]
            ),
            "reference_events": reference,
            "predicted_events_at_frozen_threshold": predicted,
            "frozen_threshold": None,
            "event_recall": {
                f"plus_minus_{tolerance}_bp": safe_fraction(
                    int(accumulator["event_matched"][index, tolerance_index]), reference
                )
                for tolerance_index, tolerance in enumerate(TOLERANCES)
            },
            "event_precision": {
                f"plus_minus_{tolerance}_bp": safe_fraction(
                    int(accumulator["event_matched"][index, tolerance_index]), predicted
                )
                for tolerance_index, tolerance in enumerate(TOLERANCES)
            },
        }
    return {
        "examples": int(accumulator["examples"]),
        "bases": int(accumulator["bases"]),
        "region": {
            "confusion_truth_rows_prediction_columns": confusion.tolist(),
            "accuracy": safe_fraction(int(np.trace(confusion)), int(confusion.sum())),
            "macro_F1": float(np.mean(f1_values)),
            "intergenic_false_positive_rate": safe_fraction(
                int(confusion[0, 1:].sum()), int(confusion[0, :].sum())
            ),
            "classes": classes,
        },
        "boundary": boundary,
        "phase": {
            "CDS_base_accuracy": safe_fraction(accumulator["phase_correct"], accumulator["phase_total"]),
            "CDS_bases": int(accumulator["phase_total"]),
            "truth_CDS_start_accuracy": safe_fraction(
                accumulator["CDS_start_phase_correct"], accumulator["CDS_start_phase_total"]
            ),
            "truth_CDS_starts": int(accumulator["CDS_start_phase_total"]),
        },
    }


def interval_overlap(left, right):
    total = 0
    for left_start, left_end in left:
        for right_start, right_end in right:
            total += max(0, min(left_end, right_end) - max(left_start, right_start))
    return total


def truth_events(cds):
    return {
        "start": [cds[0][0]],
        "stop": [cds[-1][1] - 3],
        "donor": [left[1] for left, _right in zip(cds, cds[1:])],
        "acceptor": [right[0] - 2 for _left, right in zip(cds, cds[1:])],
    }


def event_anchor(position, event):
    if event == "stop":
        return position - 3
    if event == "acceptor":
        return position - 2
    return position


def canonical_truth(sequence, cds):
    events = truth_events(cds)
    motifs_ok = all(
        sequence[position:position + len(next(iter(m25.MOTIFS[name])))] in m25.MOTIFS[name]
        for name, positions in events.items() for position in positions
    )
    coding = "".join(sequence[start:end] for start, end in cds)
    orf_ok = (
        len(coding) >= 6
        and len(coding) % 3 == 0
        and coding.startswith("ATG")
        and coding[-3:] in m25.STOP_CODONS
        and not any(coding[offset:offset + 3] in m25.STOP_CODONS
                    for offset in range(3, len(coding) - 3, 3))
    )
    return motifs_ok and orf_ok


def build_reference_records(species, validation_references):
    records = []
    for (species_name, seqid), transcripts in validation_references.items():
        sequence = species[species_name]["seqs"][seqid]
        for transcript in transcripts:
            oriented = m25._oriented_transcript(transcript, len(sequence), transcript["strand"])
            cds = [(int(start), int(end)) for start, end, _phase in oriented["CDS"]]
            oriented_sequence = sequence if transcript["strand"] == "+" else m25.reverse_complement(sequence)
            events = truth_events(cds)
            event_positions = [position for values in events.values() for position in values]
            records.append({
                "key": (species_name, seqid, transcript["strand"], transcript["id"]),
                "species": species_name,
                "seqid": seqid,
                "strand": transcript["strand"],
                "transcript_id": transcript["id"],
                "cds": cds,
                "events": events,
                "span": (cds[0][0], cds[-1][1]),
                "canonical": canonical_truth(oriented_sequence, cds),
                "span_gt_6144": cds[-1][1] - cds[0][0] > WINDOW_BP,
                "tile_edge_within_6": any(
                    min(position % WINDOW_BP, WINDOW_BP - position % WINDOW_BP) <= RADIUS_BP
                    for position in event_positions
                ),
                "CDS_count": len(cds),
            })
    return records


def transition_upper_bounds(states, sequence, reference):
    candidates = defaultdict(list)
    for position, _left, _right, events in m25.transition_candidates(states):
        for event in events:
            candidates[event].append(event_anchor(position, event))
    transition_reachable = True
    motif_reachable = True
    for name, positions in reference["events"].items():
        anchors = candidates[name]
        transition_match = event_match_count(positions, anchors, RADIUS_BP) == len(positions)
        actual_candidates = {
            position
            for anchor in anchors
            for position in motif_positions(sequence, name, anchor)
        }
        transition_reachable &= transition_match
        motif_reachable &= transition_match and all(position in actual_candidates for position in positions)
    return bool(transition_reachable), bool(motif_reachable)


def truth_assisted_exact_chain(reference, trace, sequence, boundary_logits, phase_logits, thresholds):
    if trace is None or len(trace["runs"]) != len(reference["cds"]):
        return False
    for name, truth_positions in reference["events"].items():
        candidate_rows = trace["candidate_positions"][name]
        if len(candidate_rows) != len(truth_positions):
            return False
        if any(truth not in candidates for truth, candidates in zip(truth_positions, candidate_rows)):
            return False

    boundary_probability = m25._sigmoid(boundary_logits)
    for index, name in enumerate(m25.BOUNDARY_NAMES):
        if any(boundary_probability[position, index] < float(thresholds[name])
               for position in reference["events"][name]):
            return False

    phase_class = np.asarray(phase_logits).argmax(axis=-1)
    coding_length = 0
    for index, (start, end) in enumerate(reference["cds"]):
        expected = 0 if index == 0 else (3 - coding_length % 3) % 3
        if int(phase_class[start]) != expected + 1:
            return False
        coding_length += end - start
    return canonical_truth(sequence, reference["cds"])


def assign_references(references, states_by_key, traces_by_key, sequences_by_key,
                      scores_by_key, thresholds):
    assignments = {}
    matched_lineages = {}
    upper = Counter()
    overlap_multiplicity = {}

    grouped = defaultdict(list)
    for reference in references:
        grouped[(reference["species"], reference["seqid"], reference["strand"])].append(reference)

    for key, group in grouped.items():
        states = states_by_key[key]
        traces = traces_by_key[key]
        sequence = sequences_by_key[key]
        for reference in group:
            transition, motif = transition_upper_bounds(states, sequence, reference)
            reference["transition_reachable"] = transition
            reference["motif_reachable"] = motif
            upper["transition_reachable"] += transition
            upper["motif_reachable"] += motif

        pairs = []
        for reference in group:
            covering = []
            for trace in traces:
                block_start, block_end = trace["block_span"]
                if block_start <= reference["span"][0] and block_end >= reference["span"][1]:
                    overlap = interval_overlap(reference["cds"], trace["runs"])
                    covering.append(trace["lineage_id"])
                    pairs.append((-overlap, -(block_end - block_start), reference["transcript_id"],
                                  trace["lineage_id"], reference, trace))
            overlap_multiplicity[reference["key"]] = len(covering)

        used_reference = set()
        used_lineage = set()
        for _negative_overlap, _negative_span, _tx_id, lineage_id, reference, trace in sorted(pairs):
            if reference["key"] in used_reference or lineage_id in used_lineage:
                continue
            used_reference.add(reference["key"])
            used_lineage.add(lineage_id)
            matched_lineages[reference["key"]] = trace

        for reference in group:
            trace = matched_lineages.get(reference["key"])
            boundary_logits, phase_logits = scores_by_key[key]
            truth_assisted = truth_assisted_exact_chain(
                reference, trace, sequence, boundary_logits, phase_logits, thresholds
            )
            reference["truth_assisted_exact_chain"] = truth_assisted
            upper["truth_assisted_exact_chain"] += truth_assisted

            truth_cds_overlap = sum(
                int((states[start:end] == m25.C).sum()) for start, end in reference["cds"]
            )
            if truth_cds_overlap == 0:
                stage = "region_state_path"
            elif trace is None:
                stage = "non_intergenic_block"
            elif len(trace["runs"]) != len(reference["cds"]):
                stage = "ordered_CDS_runs"
            elif trace["failure_stage"] == "legal_terminal_transitions":
                stage = "legal_terminal_transitions"
            elif trace["failure_stage"] == "start_stop_motif_candidates":
                stage = "start_stop_motif_candidates"
            elif trace["failure_stage"] == "ordered_donor_acceptor_candidates":
                stage = "ordered_donor_acceptor_candidates"
            elif any(
                truth not in candidates
                for name in ("start", "stop")
                for truth, candidates in zip(reference["events"][name], trace["candidate_positions"][name])
            ):
                stage = "start_stop_motif_candidates"
            elif any(
                truth not in candidates
                for name in ("donor", "acceptor")
                for truth, candidates in zip(reference["events"][name], trace["candidate_positions"][name])
            ):
                stage = "ordered_donor_acceptor_candidates"
            elif any(trace["chosen_positions"][name] != reference["events"][name]
                     for name in m25.BOUNDARY_NAMES):
                stage = "learned_boundary_choice"
            elif trace["failure_stage"] == "phase_check":
                stage = "phase_check"
            elif trace["failure_stage"] == "complete_ORF_internal_stop_check":
                stage = "complete_ORF_internal_stop_check"
            elif trace["failure_stage"] == "boundary_threshold_filter":
                stage = "boundary_threshold_filter"
            elif trace.get("emitted_model") and tuple(trace["emitted_model"]["cds"]) == tuple(reference["cds"]):
                stage = "emitted_exact_chain"
            else:
                stage = "learned_boundary_choice"
            assignments[reference["key"]] = stage

    if len(assignments) != len(references):
        raise AssertionError("reference assignment is incomplete")
    return assignments, upper, overlap_multiplicity, matched_lineages


def count_by_strata(references, assignments):
    strata = defaultdict(Counter)
    for reference in references:
        labels = [
            f"species:{reference['species']}",
            "single_CDS" if reference["CDS_count"] == 1 else "multi_CDS",
            ("CDS_count:1" if reference["CDS_count"] == 1 else
             "CDS_count:2" if reference["CDS_count"] == 2 else
             "CDS_count:3-5" if reference["CDS_count"] <= 5 else "CDS_count:6+"),
            "span_gt_6144" if reference["span_gt_6144"] else "span_lte_6144",
            "tile_edge_within_6" if reference["tile_edge_within_6"] else "not_near_tile_edge",
            "canonical" if reference["canonical"] else "noncanonical",
        ]
        for label in labels:
            strata[label][assignments[reference["key"]]] += 1
            strata[label]["reference_total"] += 1
    return {label: dict(counts) for label, counts in sorted(strata.items())}


def model_transcripts(predictions):
    transcripts = []
    serial = 0
    for (species_name, seqid), models in predictions.items():
        for model in models:
            serial += 1
            transcripts.append({
                "id": f"prediction_{serial}",
                "gene_id": f"prediction_{serial}",
                "species": species_name,
                "seqid": seqid,
                "strand": model["strand"],
                "CDS": [(start, end, str(phase))
                        for (start, end), phase in zip(model["cds"], model["phase"])],
                "exon": [],
                "partial": False,
            })
    return transcripts


def prediction_error_summary(references, prediction_transcripts):
    reference_transcripts = [tx for rows in references.values() for tx in rows]
    ref_exact = {(tx["seqid"], tx["strand"], tuple((a, b) for a, b, _p in tx["CDS"])) for tx in reference_transcripts}
    pred_exact = {(tx["seqid"], tx["strand"], tuple((a, b) for a, b, _p in tx["CDS"])) for tx in prediction_transcripts}
    ref_no_strand = {(seqid, chain) for seqid, _strand, chain in ref_exact}
    wrong_strand = sum((seqid, chain) in ref_no_strand and (seqid, strand, chain) not in ref_exact
                       for seqid, strand, chain in pred_exact)
    return {
        "reference_chains": len(ref_exact),
        "predicted_chains": len(pred_exact),
        "exact_chains": len(ref_exact & pred_exact),
        "missing_chains": len(ref_exact - pred_exact),
        "false_positive_chains": len(pred_exact - ref_exact),
        "exact_coordinates_wrong_strand": int(wrong_strand),
        "overlap_split_fusion": m25_eval.m24.overlap_degrees(prediction_transcripts, reference_transcripts),
        "matched_gene_strand": m25_eval.matched_gene_strand_accuracy(reference_transcripts, prediction_transcripts),
        "exact_matched_CDS_phase": m25_eval.phase_accuracy(reference_transcripts, prediction_transcripts),
        "boundary_offsets": m25_eval.boundary_diagnostics(reference_transcripts, prediction_transcripts),
    }


def audit_gff3(path, sequences, lengths, expected_models):
    genes = {}
    transcripts = {}
    children = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_no}: expected 9 columns")
            seqid, _source, feature, start_text, end_text, _score, strand, phase, attrs_text = fields
            start, end = int(start_text) - 1, int(end_text)
            attrs = m25_eval.m24.parse_attrs(attrs_text)
            if seqid not in lengths or start < 0 or end <= start or end > lengths[seqid] or strand not in {"+", "-"}:
                raise ValueError(f"{path}:{line_no}: invalid emitted coordinate or strand")
            if feature == "gene":
                genes[attrs["ID"]] = (seqid, start, end, strand)
            elif feature == "mRNA":
                transcripts[attrs["ID"]] = (attrs["Parent"], seqid, start, end, strand)
            elif feature in {"CDS", "start_codon", "stop_codon"}:
                children[attrs["Parent"]].append((feature, seqid, start, end, strand, phase))

    integrity_failures = []
    for transcript_id, (gene_id, seqid, start, end, strand) in transcripts.items():
        if gene_id not in genes or genes[gene_id] != (seqid, start, end, strand):
            integrity_failures.append({"transcript": transcript_id, "reason": "parent_linkage"})
        if not children[transcript_id]:
            integrity_failures.append({"transcript": transcript_id, "reason": "missing_children"})
    for transcript_id in children:
        if transcript_id not in transcripts:
            integrity_failures.append({"transcript": transcript_id, "reason": "orphan_child"})
    if integrity_failures:
        raise AssertionError(f"emitted GFF3 linkage audit failed: {integrity_failures}")

    annotation = parse_annotation(str(path), lengths, protein_coding_only=False)
    parsed = primary_transcripts(annotation)
    n_checked = len(parsed)
    if n_checked != expected_models or len(transcripts) != expected_models or len(genes) != expected_models:
        raise AssertionError("emitted GFF3 model count does not reconcile")

    codons = m25_eval.codon_features(path)
    transcript_ledger = []
    component_failure_counts = Counter()
    failures = []
    for transcript in parsed:
        components = {
            "parent_linkage": True,
            **m25_eval.structural_validity_components(transcript, sequences, codons),
        }
        failed_components = [
            name for name in VALIDITY_COMPONENTS if not components[name]
        ]
        for name in failed_components:
            component_failure_counts[name] += 1
        transcript_ledger.append({
            "transcript": transcript["id"],
            "valid": not failed_components,
            "components": components,
            "failed_components": failed_components,
        })
        if failed_components:
            failures.append({
                "transcript": transcript["id"],
                "failed_components": failed_components,
            })

    valid_transcripts = sum(row["valid"] for row in transcript_ledger)
    return {
        "n_emitted": expected_models,
        "n_checked": n_checked,
        "audit_coverage": safe_fraction(n_checked, expected_models) if expected_models else 1.0,
        "complete_empty_audit": expected_models == 0 and n_checked == 0,
        "valid_transcripts": valid_transcripts,
        "invalid_transcripts": n_checked - valid_transcripts,
        "validity_fraction": (safe_fraction(valid_transcripts, n_checked)
                              if n_checked else "not_applicable"),
        "component_failure_counts": {
            name: int(component_failure_counts[name]) for name in VALIDITY_COMPONENTS
        },
        "failures": failures,
        "transcript_ledger": transcript_ledger,
    }


def load_species(root, config):
    primary_seqids = config["data"]["primary_chromosome_seqids"]
    species = {}
    for relative_path in config["data"]["development_species"]:
        species_path = root / relative_path
        name = species_path.name
        if name not in {"arabidopsis_thaliana", "oryza_sativa"}:
            raise ValueError(f"Stage 1 development species is not frozen: {name}")
        all_sequences = screen_data.read_fasta(str(species_path / "genome.fa"))
        all_splits = screen_data.assign_splits(list(all_sequences))
        sequences = m25.select_primary_chromosomes(all_sequences, primary_seqids[name])
        lengths = {seqid: len(sequence) for seqid, sequence in sequences.items()}
        annotation = parse_annotation(str(species_path / "reference.gff3"), lengths, protein_coding_only=True)
        transcripts = m25._primary_with_partial(annotation)
        by_seqid = defaultdict(list)
        for transcript in transcripts:
            by_seqid[transcript["seqid"]].append(transcript)
        species[name] = {
            "seqs": sequences,
            "lengths": lengths,
            "annotation": annotation,
            "transcripts": transcripts,
            "transcripts_by_seqid": by_seqid,
            "splits": {seqid: all_splits[seqid] for seqid in sequences},
        }
    return species


def validation_truth(species):
    references = {}
    lengths = {}
    for species_name, record in species.items():
        by_seqid = defaultdict(list)
        for transcript in primary_transcripts(record["annotation"]):
            if record["splits"][transcript["seqid"]] == "val":
                by_seqid[transcript["seqid"]].append(transcript)
        for seqid, split in record["splits"].items():
            if split == "val":
                references[(species_name, seqid)] = by_seqid[seqid]
                lengths[(species_name, seqid)] = record["lengths"][seqid]
    return references, lengths


def load_inference_model(config):
    import torch
    from peft import LoraConfig, get_peft_model
    from transformers import AutoModelForTokenClassification, AutoTokenizer

    model_name = config["model"]["model_name"]
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    full_model = AutoModelForTokenClassification.from_pretrained(
        model_name, trust_remote_code=True, attn_implementation="sdpa"
    )
    k = int(getattr(full_model, "k", getattr(full_model.config, "k", 6)))
    if k != 6:
        raise ValueError(f"frozen model reports k={k}, expected 6")
    backbone = full_model.model
    del full_model
    backbone.gradient_checkpointing_enable()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device != "cuda":
        raise RuntimeError("Stage 1 requires the frozen GPU inference path")
    backbone.to(device=device, dtype=torch.bfloat16)
    for parameter in backbone.parameters():
        parameter.requires_grad_(False)
    lora = config["lora"]
    backbone = get_peft_model(backbone, LoraConfig(
        r=int(lora["r"]),
        lora_alpha=int(lora["alpha"]),
        target_modules=list(lora["targets"]),
        lora_dropout=float(lora["dropout"]),
        bias="none",
    ))

    class StructuralHeads(torch.nn.Module):
        def __init__(self, hidden_size):
            super().__init__()
            channels = hidden_size + 5
            self.region = _ConvLSTMHead(channels, n_classes=3)
            self.boundary = torch.nn.Conv1d(channels, 4, kernel_size=9, padding=4)
            self.phase = torch.nn.Conv1d(channels, 4, kernel_size=1)

        def forward(self, per_base, nucleotide):
            features = torch.cat((per_base, nucleotide.to(per_base.dtype)), dim=-1).transpose(1, 2)
            return (self.region(features), self.boundary(features).transpose(1, 2),
                    self.phase(features).transpose(1, 2))

    heads = StructuralHeads(int(backbone.config.hidden_size)).to(device)

    def forward(ids, attention, nucleotide):
        hidden = backbone(input_ids=ids, attention_mask=attention).last_hidden_state
        per_base = hidden.repeat_interleave(k, dim=1)
        return heads(per_base, nucleotide)

    return tokenizer, backbone, heads, forward, device, k


def load_checkpoint(path, backbone, heads):
    import torch

    checkpoint = torch.load(path, map_location="cpu")
    heads.load_state_dict(checkpoint["heads"])
    current = backbone.state_dict()
    current.update(checkpoint["lora"])
    backbone.load_state_dict(current)
    backbone.eval()
    heads.eval()


def predict_sequence(sequence, tokenizer, forward, device, k):
    import torch

    score = [np.zeros((len(sequence), width), dtype=np.float16) for width in (3, 4, 4)]
    with torch.no_grad():
        for start in range(0, len(sequence), WINDOW_BP):
            real_end = min(start + WINDOW_BP, len(sequence))
            padded = sequence[start:real_end] + "A" * (WINDOW_BP - (real_end - start))
            ids = _tokenize_window(tokenizer, _clean(padded), WINDOW_BP, k).unsqueeze(0).to(device)
            attention = torch.ones_like(ids)
            nucleotide = torch.from_numpy(m25._one_hot(padded)).unsqueeze(0).to(device)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                outputs = forward(ids, attention, nucleotide)
            width = real_end - start
            for target, output in zip(score, outputs):
                target[start:real_end] = output[0, :width].float().cpu().numpy().astype(np.float16)
    return score


def predict_train_raw(dataset, thresholds, forward, device):
    import torch

    accumulator = new_raw_accumulator()
    with torch.no_grad():
        for index in range(len(dataset)):
            ids, attention, nucleotide, region, boundary, phase, structural_mask = dataset[index]
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                outputs = forward(
                    ids.unsqueeze(0).to(device), attention.unsqueeze(0).to(device),
                    nucleotide.unsqueeze(0).to(device),
                )
            arrays = [output[0].float().cpu().numpy().astype(np.float16) for output in outputs]
            add_raw_metrics(accumulator, *arrays, {
                "region": region.numpy(),
                "boundary": boundary.numpy(),
                "phase": phase.numpy(),
                "structural_mask": structural_mask.numpy(),
            }, thresholds)
            if (index + 1) % 128 == 0:
                print(f"train_raw_examples={index + 1}", flush=True)
    result = finalize_raw_metrics(accumulator)
    for name in m25.BOUNDARY_NAMES:
        result["boundary"][name]["frozen_threshold"] = float(thresholds[name])
    return result


def run_epoch(epoch, row, checkpoint_path, species, validation_references, validation_lengths,
              reference_records, train_dataset, tokenizer, backbone, heads, forward, device, k, out_dir):
    thresholds = {name: float(row[name]) for name in ("region",) + m25.BOUNDARY_NAMES}
    load_checkpoint(checkpoint_path, backbone, heads)
    print(f"epoch={epoch} train raw-head inference", flush=True)
    train_raw = predict_train_raw(train_dataset, thresholds, forward, device)

    validation_raw_accumulator = new_raw_accumulator()
    validation_scores = {}
    print(f"epoch={epoch} validation chromosome inference", flush=True)
    for species_name, record in species.items():
        for seqid, split in record["splits"].items():
            if split != "val":
                continue
            sequence = record["seqs"][seqid]
            validation_scores[(species_name, seqid)] = {}
            for strand, oriented in (("+", sequence), ("-", m25.reverse_complement(sequence))):
                scores = predict_sequence(oriented, tokenizer, forward, device, k)
                validation_scores[(species_name, seqid)][strand] = scores
                targets = m25.build_strand_targets(len(sequence), record["transcripts_by_seqid"][seqid], strand)
                add_raw_metrics(validation_raw_accumulator, *scores, targets, thresholds)
            print(f"epoch={epoch} validation_done={species_name}/{seqid}", flush=True)
    validation_raw = finalize_raw_metrics(validation_raw_accumulator)
    for name in m25.BOUNDARY_NAMES:
        validation_raw["boundary"][name]["frozen_threshold"] = float(thresholds[name])

    predictions = {}
    states_by_key = {}
    traces_by_key = {}
    sequences_by_key = {}
    scores_by_key = {}
    lineage_counts = Counter()
    phase_checks = []
    epoch_dir = out_dir / f"epoch_{epoch}"
    epoch_dir.mkdir(parents=True, exist_ok=True)
    lineage_path = epoch_dir / "candidate_lineages.jsonl"
    with lineage_path.open("w", encoding="utf-8") as lineage_handle:
        for key, orientation_scores in validation_scores.items():
            species_name, seqid = key
            sequence = species[species_name]["seqs"][seqid]
            models = []
            for strand, oriented in (("+", sequence), ("-", m25.reverse_complement(sequence))):
                region_score, boundary_score, phase_score = orientation_scores[strand]
                states, traces, prefilter, emitted = trace_decode_orientation(
                    oriented, region_score, boundary_score, phase_score, thresholds,
                    reverse_mapped=strand == "-",
                )
                trace_key = (species_name, seqid, strand)
                states_by_key[trace_key] = states
                traces_by_key[trace_key] = traces
                sequences_by_key[trace_key] = oriented
                scores_by_key[trace_key] = (boundary_score, phase_score)
                lineage_counts["non_intergenic_blocks"] += len(traces)
                lineage_counts["prefilter_models"] += len(prefilter)
                lineage_counts["emitted_models"] += len(emitted)
                for trace in traces:
                    lineage_counts[trace["failure_stage"] or "emitted"] += 1
                    phase_checks.extend(trace["phase_checks"])
                    payload = {
                        "species": species_name,
                        "seqid": seqid,
                        "strand": strand,
                        **trace,
                    }
                    lineage_handle.write(json.dumps(payload, sort_keys=True) + "\n")
                models.extend(m25._map_model_from_orientation(model, len(sequence), strand) for model in emitted)
            predictions[key] = models

    metrics = m25._validation_metrics(predictions, validation_references, validation_lengths)
    frozen_metrics = row["metrics"]
    reproduction = {
        name: {
            "frozen": float(frozen_metrics[name]),
            "replayed": float(metrics[name]),
            "absolute_error": abs(float(metrics[name]) - float(frozen_metrics[name])),
        }
        for name in ("exact_CDS_interval_F1", "exact_CDS_chain_F1", "intergenic_FPR", "gene_count_ratio")
    }
    if any(item["absolute_error"] > 1e-5 for item in reproduction.values()):
        raise AssertionError(f"epoch {epoch} frozen aggregate reproduction failed: {reproduction}")

    assignments, upper, overlap_multiplicity, matched_lineages = assign_references(
        reference_records, states_by_key, traces_by_key, sequences_by_key, scores_by_key, thresholds
    )
    stage_counts = Counter(assignments.values())
    if sum(stage_counts.values()) != EXPECTED_REFERENCE_CHAINS:
        raise AssertionError("reference attrition counts do not sum to 6,450")

    attrition_path = epoch_dir / "reference_attrition.tsv"
    with attrition_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "species", "seqid", "strand", "transcript_id", "earliest_stage", "CDS_count",
            "canonical", "span_gt_6144", "tile_edge_within_6", "covering_lineages",
            "transition_reachable", "motif_reachable", "truth_assisted_exact_chain",
        ])
        for reference in sorted(reference_records, key=lambda item: item["key"]):
            writer.writerow([
                reference["species"], reference["seqid"], reference["strand"], reference["transcript_id"],
                assignments[reference["key"]], reference["CDS_count"], int(reference["canonical"]),
                int(reference["span_gt_6144"]), int(reference["tile_edge_within_6"]),
                overlap_multiplicity[reference["key"]], int(reference["transition_reachable"]),
                int(reference["motif_reachable"]), int(reference["truth_assisted_exact_chain"]),
            ])

    models_by_seqid = defaultdict(list)
    for (_species_name, seqid), models in predictions.items():
        models_by_seqid[seqid].extend(models)
    gff_path = epoch_dir / "replayed_predictions.gff3"
    m25._write_gff3(str(gff_path), models_by_seqid, "M25R_development_replay")
    flat_predictions = m25.predictions_flat(predictions)
    sequences = {seqid: record["seqs"][seqid] for record in species.values() for seqid in record["seqs"]}
    lengths = {seqid: len(sequence) for seqid, sequence in sequences.items()}
    validity = audit_gff3(gff_path, sequences, lengths, len(flat_predictions))
    transcript_ledger = validity.pop("transcript_ledger")
    if len(transcript_ledger) != validity["n_checked"]:
        raise AssertionError("structural validity ledger does not reconcile")
    validity_path = epoch_dir / "structural_validity.jsonl"
    with validity_path.open("w", encoding="utf-8") as handle:
        for row in transcript_ledger:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    if validity["audit_coverage"] != 1.0:
        raise AssertionError(f"independent GFF3 validity audit failed: {validity}")

    terminal_status_counts = {
        key: int(value)
        for key, value in lineage_counts.items()
        if key not in {"non_intergenic_blocks", "prefilter_models", "emitted_models"}
    }
    if sum(terminal_status_counts.values()) != lineage_counts["non_intergenic_blocks"]:
        raise AssertionError("candidate lineages do not have exactly one terminal status")
    if lineage_counts["emitted_models"] != len(flat_predictions):
        raise AssertionError("emitted lineage count does not match final GFF3 model count")

    prediction_transcripts = model_transcripts(predictions)
    errors = prediction_error_summary(validation_references, prediction_transcripts)
    candidate_phase_accuracy = safe_fraction(
        sum(item["predicted"] == item["expected"] for item in phase_checks), len(phase_checks)
    )
    examples = defaultdict(list)
    for reference in sorted(reference_records, key=lambda item: item["key"]):
        stage = assignments[reference["key"]]
        if len(examples[stage]) < 5:
            examples[stage].append({
                "species": reference["species"], "seqid": reference["seqid"],
                "strand": reference["strand"], "transcript_id": reference["transcript_id"],
                "CDS_count": reference["CDS_count"], "canonical": reference["canonical"],
            })

    canonical_references = [reference for reference in reference_records if reference["canonical"]]
    per_species_metrics = {}
    for species_name in sorted(species):
        species_predictions = {key: value for key, value in predictions.items() if key[0] == species_name}
        species_references = {key: value for key, value in validation_references.items() if key[0] == species_name}
        species_lengths = {key: value for key, value in validation_lengths.items() if key[0] == species_name}
        per_species_metrics[species_name] = m25._validation_metrics(
            species_predictions, species_references, species_lengths
        )

    result = {
        "epoch": epoch,
        "checkpoint": str(checkpoint_path),
        "frozen_tuple": {
            name: row[name] for name in ("epoch", "region", "start", "stop", "donor", "acceptor", "enumeration_order")
        },
        "reproduction": reproduction,
        "train_raw_heads": train_raw,
        "validation_raw_heads": validation_raw,
        "reference_attrition": {
            "reference_total": len(reference_records),
            "stage_counts": {stage: int(stage_counts.get(stage, 0)) for stage in STAGES},
            "stage_fractions": {stage: stage_counts.get(stage, 0) / len(reference_records) for stage in STAGES},
            "canonical_reference_count": sum(reference["canonical"] for reference in reference_records),
            "unsupported_by_frozen_canonical_decoder": sum(not reference["canonical"] for reference in reference_records),
            "strata": count_by_strata(reference_records, assignments),
            "representative_examples": dict(examples),
        },
        "post_hoc_upper_bounds": {
            name: {
                "chains": int(upper[name]),
                "fraction_of_R_all": upper[name] / len(reference_records),
                "canonical_chains": sum(reference[name] for reference in canonical_references),
                "fraction_of_R_canonical": (sum(reference[name] for reference in canonical_references)
                                            / len(canonical_references) if canonical_references else "not_applicable"),
            }
            for name in ("transition_reachable", "motif_reachable", "truth_assisted_exact_chain")
        },
        "candidate_lineage_reconciliation": {
            "non_intergenic_blocks": int(lineage_counts["non_intergenic_blocks"]),
            "terminal_status_counts": terminal_status_counts,
            "prefilter_models": int(lineage_counts["prefilter_models"]),
            "boundary_filter_rejections": int(lineage_counts["boundary_threshold_filter"]),
            "emitted_lineages": int(lineage_counts["emitted_models"]),
            "final_GFF3_models": len(flat_predictions),
            "unreconciled_lineages": 0,
        },
        "candidate_start_expected_phase_accuracy": candidate_phase_accuracy,
        "per_species_metrics": per_species_metrics,
        "prediction_errors": errors,
        "independent_GFF3_validity": validity,
        "paths": {
            "reference_attrition": str(attrition_path),
            "candidate_lineages": str(lineage_path),
            "replayed_predictions": str(gff_path),
            "structural_validity": str(validity_path),
        },
    }
    save_json_atomic(epoch_dir / "diagnostic.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--experiment-id", default=EXPERIMENT_ID)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--m25r-output", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    import torch
    import yaml

    torch.manual_seed(0)
    np.random.seed(0)
    root = args.root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    with args.config.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if config["seed"] != 0 or config["training"]["epochs"] != 3 or config["model"]["window_bp"] != WINDOW_BP:
        raise ValueError("M25R frozen seed/epoch/window contract changed")

    grid_path = args.m25r_output / "validation_grid_diagnostics.json"
    with grid_path.open(encoding="utf-8") as handle:
        grid = json.load(handle)
    rows = select_epoch_rows(grid)
    species = load_species(root, config)
    validation_references, validation_lengths = validation_truth(species)
    reference_records = build_reference_records(species, validation_references)
    reference_chain_count = len({
        (species_name, seqid, tx["strand"], tuple((a, b) for a, b, _phase in tx["CDS"]))
        for (species_name, seqid), transcripts in validation_references.items() for tx in transcripts
    })
    if len(reference_records) != EXPECTED_REFERENCE_CHAINS or reference_chain_count != EXPECTED_REFERENCE_CHAINS:
        raise ValueError(
            f"frozen reference universe changed: records={len(reference_records)} chains={reference_chain_count}"
        )

    tokenizer, backbone, heads, forward, device, k = load_inference_model(config)
    train_dataset = m25.OrientationWindowDataset(
        species, config["data"]["train_split"], tokenizer, WINDOW_BP,
        float(config["training"]["sample_fraction"]), int(config["seed"]),
        int(config["training"]["train_window_cap"]),
    )
    if len(train_dataset) != EXPECTED_TRAIN_WINDOWS:
        raise ValueError(f"frozen training example set changed: {len(train_dataset)}")

    inputs = {
        "experiment_id": args.experiment_id,
        "source_experiment": config["exp_id"],
        "config": str(args.config.resolve()),
        "validation_grid": str(grid_path.resolve()),
        "development_species": [str((root / path).resolve()) for path in config["data"]["development_species"]],
        "reference_chains": len(reference_records),
        "train_examples": len(train_dataset),
        "epochs": {
            str(epoch): {
                "checkpoint": str((args.m25r_output / "checkpoints" / f"epoch_{epoch}.pt").resolve()),
                "grid_row_enumeration_order": int(rows[epoch]["enumeration_order"]),
            }
            for epoch in (1, 2, 3)
        },
        "setaria_files_read": False,
        "weights_updated": False,
        "threshold_or_decoder_search": False,
    }
    save_json_atomic(out_dir / "resolved_inputs.json", inputs)

    epoch_results = []
    for epoch in (1, 2, 3):
        checkpoint_path = (args.m25r_output / "checkpoints" / f"epoch_{epoch}.pt").resolve()
        epoch_result = run_epoch(
            epoch, rows[epoch], checkpoint_path, species, validation_references, validation_lengths,
            reference_records, train_dataset, tokenizer, backbone, heads, forward, device, k, out_dir,
        )
        if not (out_dir / f"epoch_{epoch}" / "diagnostic.json").is_file():
            raise AssertionError(f"epoch {epoch} diagnostic was not persisted")
        epoch_results.append(epoch_result)

    stage_integrity = {
        "aggregate_reproduction_within_1e-5": all(
            item["absolute_error"] <= 1e-5
            for result in epoch_results for item in result["reproduction"].values()
        ),
        "all_6450_references_accounted_each_epoch": all(
            sum(result["reference_attrition"]["stage_counts"].values()) == EXPECTED_REFERENCE_CHAINS
            for result in epoch_results
        ),
        "validity_audit_coverage_100_percent": all(
            result["independent_GFF3_validity"]["audit_coverage"] == 1.0 for result in epoch_results
        ),
        "setaria_files_read": False,
        "weights_updated": False,
        "threshold_or_decoder_search": False,
    }
    any_invalid = any(
        result["independent_GFF3_validity"]["invalid_transcripts"] > 0
        for result in epoch_results
    )
    result = {
        "experiment_id": args.experiment_id,
        "status": "COMPLETED_STAGE1_REVIEW_REQUIRED",
        "scientific_status": ("SCIENTIFIC_NO_GO_INVALID_STRUCTURES" if any_invalid
                              else "DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED"),
        "stage_integrity": stage_integrity,
        "epochs": epoch_results,
        "scientific_next_action": "STOP_FOR_REVIEW",
    }
    save_json_atomic(out_dir / "stage1_diagnostic.json", result)
    (out_dir / "STATUS").write_text("COMPLETED_STAGE1_REVIEW_REQUIRED\n", encoding="utf-8")
    print("COMPLETED_STAGE1_REVIEW_REQUIRED", flush=True)


if __name__ == "__main__":
    main()
