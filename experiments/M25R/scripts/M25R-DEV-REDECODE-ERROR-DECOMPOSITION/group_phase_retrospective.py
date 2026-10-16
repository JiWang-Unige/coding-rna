#!/usr/bin/env python3
"""Retrospective grouping/phase evidence from saved R4 traces; no model inference."""
import csv
import json
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import posthoc_r4_case_correction as case

D = case.diag
ROOT = Path(__file__).resolve().parents[3]
OUT_ID = "M25R-R4-GROUP-PHASE-RETROSPECTIVE"
NAMES = D.m25.BOUNDARY_NAMES


def chain(parts):
    return tuple(tuple(p) for p in parts)


def chosen_cds(trace):
    p = trace["chosen_positions"]
    n = len(trace["runs"])
    if n == 0 or len(p["start"]) != 1 or len(p["stop"]) != 1 or any(
            len(p[name]) != n - 1 for name in ("donor", "acceptor")):
        return ()
    starts = p["start"] + [a + 2 for a in p["acceptor"]]
    ends = p["donor"] + [p["stop"][0] + 3]
    if any(a >= b for a, b in zip(starts, ends)):
        return ()
    return tuple(zip(starts, ends))


def contains_truth(ref, trace):
    return len(trace["runs"]) == len(ref["cds"]) and all(
        len(trace["candidate_positions"][name]) == len(ref["events"][name])
        and all(p in options for p, options in zip(ref["events"][name], trace["candidate_positions"][name]))
        for name in NAMES)


def original_matching(refs, traces):
    """Reproduce R4 exact-emitted-first, then frozen greedy covering-pair order."""
    ordered = sorted(traces, key=lambda t: t["block_span"])
    starts = [t["block_span"][0] for t in ordered]
    emitted = defaultdict(list)
    for t in traces:
        if t.get("emitted_model"):
            emitted[chain(t["emitted_model"]["cds"])].append(t)
    exact, pairs, cover = [], [], {}
    for r in refs:
        i = bisect_right(starts, r["span"][0]) - 1
        t = ordered[i] if i >= 0 and ordered[i]["block_span"][1] >= r["span"][1] else None
        cover[r["key"]] = t
        if t is not None:
            a, b = t["block_span"]
            pairs.append((-D.interval_overlap(r["cds"], t["runs"]), -(b - a),
                          r["transcript_id"], t["lineage_id"], r, t))
        for t in emitted[chain(r["cds"])]:
            exact.append((r["transcript_id"], t["lineage_id"], r, t))
    matched, used = {}, set()
    for _id, lid, r, t in sorted(exact, key=lambda x: (x[0], x[1])):
        if r["key"] not in matched and lid not in used:
            matched[r["key"]] = t
            used.add(lid)
    for _overlap, _span, _id, lid, r, t in sorted(pairs, key=lambda x: x[:4]):
        if r["key"] not in matched and lid not in used:
            matched[r["key"]] = t
            used.add(lid)
    return matched, cover


def original_stage(ref, trace, cds_support):
    if not cds_support:
        return "region_state_path"
    if trace is None:
        return "non_intergenic_block"
    if len(trace["runs"]) != len(ref["cds"]):
        return "ordered_CDS_runs"
    fail = trace["failure_stage"]
    if fail in {"legal_terminal_transitions", "start_stop_motif_candidates", "ordered_donor_acceptor_candidates"}:
        return fail
    for names, stage in [(("start", "stop"), "start_stop_motif_candidates"),
                         (("donor", "acceptor"), "ordered_donor_acceptor_candidates")]:
        if any(p not in candidates for name in names
               for p, candidates in zip(ref["events"][name], trace["candidate_positions"][name])):
            return stage
    if any(trace["chosen_positions"][name] != ref["events"][name] for name in NAMES):
        return "learned_boundary_choice"
    if fail in {"phase_check", "complete_ORF_internal_stop_check", "boundary_threshold_filter"}:
        return fail
    if trace.get("emitted_model") and chain(trace["emitted_model"]["cds"]) == chain(ref["cds"]):
        return "emitted_exact_chain"
    return "learned_boundary_choice"


def proper_subruns(ref, trace, sequence):
    """Candidate-coordinate feasibility only, not a re-decode or score/phase oracle."""
    runs = trace["runs"]
    n = len(ref["cds"])
    if len(runs) <= n:
        return False
    for offset in range(len(runs) - n + 1):
        subset = runs[offset:offset + n]
        if subset[0][0] == 0 or subset[-1][1] == len(sequence):
            continue
        anchors = D.truth_events(subset)
        if all(len(anchors[name]) == len(ref["events"][name]) and all(
                truth in D.motif_positions(sequence, name, anchor)
                for truth, anchor in zip(ref["events"][name], anchors[name])) for name in NAMES):
            return True
    return False


def validate_phase_trace(trace):
    cds = chosen_cds(trace)
    checks = trace["phase_checks"]
    if not cds or not checks or len(checks) > len(cds):
        raise AssertionError("phase-rejected trace lacks complete chosen CDS or checked prefix")
    cumulative = 0
    for i, check in enumerate(checks):
        start, end = cds[i]
        if check["position"] != start or check["expected"] != 1 + ((3 - cumulative % 3) % 3):
            raise AssertionError("saved phase expectation inconsistent with chosen CDS prefix")
        if (check["predicted"] != check["expected"]) != (i == len(checks) - 1):
            raise AssertionError("phase checks do not stop at first rejection")
        cumulative += end - start
    return checks[-1], len(checks) - 1


def phase_evidence(trace, ref):
    failure, index = validate_phase_trace(trace)
    if ref is None:
        label = "unassigned_reference_unidentifiable"
        truth = None
    else:
        truth = ref["phase_by_start"].get(failure["position"])
        exact = chosen_cds(trace) == chain(ref["cds"])
        if truth is None:
            label = "off_reference_CDS_start_unidentifiable_phase"
        elif truth != failure["expected"]:
            if exact:
                label = "exact_chain_annotation_decoder_phase_disagreement"
            elif truth == failure["predicted"]:
                label = "wrong_chain_frame_shift_head_agrees_reference"
            else:
                label = "wrong_chain_frame_shift_and_head_disagreement"
        elif exact:
            label = "exact_chain_correct_start_phase_disagreement"
        else:
            label = "wrong_chain_correct_local_start_phase_disagreement"
    return {"label": label, "failed_CDS_index": index, "failure": failure, "reference_phase_code": truth}


def token_phase_conflict(ref, position, sequence):
    """A same-token, same-base pair has identical input to the frozen 1x1 phase head."""
    for (a, b), code in zip(ref["cds"], ref["phase_codes"]):
        if a <= position < b and code is not None:
            target = 1 + ((position - a + code - 1) % 3)
            left, right = max(a, position // 6 * 6), min(b, position // 6 * 6 + 6)
            return any(sequence[q] == sequence[position] and
                       1 + ((q - a + code - 1) % 3) != target for q in range(left, right))
    return False


def primary_category(ref, emitted, exact_chosen, candidates, subruns, support, overlap_blocks, covering):
    if emitted:
        return "recovered_exact"
    if exact_chosen:
        failures = {t["failure_stage"] for t in exact_chosen}
        if len(failures) != 1:
            return "exact_candidates_mixed_terminal_failures_unidentifiable"
        fail = next(iter(failures))
        return {"phase_check": "correct_chosen_chain_phase_rejection",
                "complete_ORF_internal_stop_check": "correct_chosen_chain_ORF_rejection",
                "boundary_threshold_filter": "correct_chosen_chain_threshold_rejection"}.get(
                    fail, "correct_chosen_chain_other_terminal")
    if candidates:
        failures = {t["failure_stage"] for t in candidates}
        if failures == {"phase_check"}:
            return "true_candidate_available_wrong_choice_then_phase"
        return "true_candidate_available_wrong_choice_other_or_mixed"
    if subruns:
        return "whole_block_excludes_truth_compatible_subruns"
    if not support:
        return "no_CDS_state_support"
    if covering is None:
        if overlap_blocks >= 2:
            return "CDS_support_fragmented_across_blocks"
        return "partial_block_span_support"
    if len(covering["runs"]) > len(ref["cds"]):
        return "covering_block_extra_runs_no_compatible_subchain"
    if len(covering["runs"]) < len(ref["cds"]):
        return "covering_block_missing_runs"
    return "same_run_count_truth_boundary_candidates_unavailable"


def analyze_epoch(epoch, species, references, output):
    source = ROOT / "outputs" / case.SOURCE_ID / f"epoch_{epoch}"
    corrected = ROOT / "outputs" / case.OUTPUT_ID / f"epoch_{epoch}"
    old = case.read_json(source / "diagnostic.json")
    with (corrected / "reference_case_correction.tsv").open() as h:
        ledger = {tuple(row[k] for k in ("species", "seqid", "strand", "transcript_id")): row
                  for row in csv.DictReader(h, delimiter="\t")}
    grouped = defaultdict(list)
    for ref in references:
        grouped[ref["key"][:3]].append(ref)
    traces_by_key = defaultdict(list)
    with (source / "candidate_lineages.jsonl").open() as h:
        for line in h:
            t = json.loads(line)
            traces_by_key[(t["species"], t["seqid"], t["strand"])].append(t)
    rows, phase_rows = [], []
    for key, refs in grouped.items():
        sequence = species[key[0]]["seqs"][key[1]]
        if key[2] == "-":
            sequence = D.m25.reverse_complement(sequence)
        traces = sorted(traces_by_key[key], key=lambda t: t["block_span"])
        states = case.reconstruct_states(len(sequence), traces)
        matched, cover = original_matching(refs, traces)
        assigned = {t["lineage_id"]: r for r in refs if (t := matched.get(r["key"])) is not None}
        starts = [t["block_span"][0] for t in traces]
        ends = [t["block_span"][1] for t in traces]
        chosen_map = defaultdict(list)
        for t in traces:
            if chosen_cds(t):
                chosen_map[chosen_cds(t)].append(t)
        ref_by_chain = defaultdict(list)
        for ref in refs:
            ref_by_chain[chain(ref["cds"])].append(ref)
        for t in traces:
            if t["failure_stage"] != "phase_check":
                continue
            exact_refs = ref_by_chain.get(chosen_cds(t), [])
            ref = exact_refs[0] if len(exact_refs) == 1 else assigned.get(t["lineage_id"])
            evidence = phase_evidence(t, ref)
            phase_rows.append({"species": key[0], "seqid": key[1], "strand": key[2],
                               "lineage_id": t["lineage_id"], "reference_id": ref["transcript_id"] if ref else None,
                               "reference_relation": "exact_chosen_chain" if len(exact_refs) == 1 else "original_assignment",
                               "chosen_CDS": chosen_cds(t), **evidence})
        for ref in refs:
            saved = ledger[ref["key"]]
            if int(saved["canonical"]) != ref["canonical"]:
                raise AssertionError("corrected canonical membership changed")
            support = sum(int(np.count_nonzero(states[a:b] == D.m25.C)) for a, b in ref["cds"])
            match = matched.get(ref["key"])
            stage = original_stage(ref, match, support)
            if stage != saved["earliest_stage"]:
                raise AssertionError(f"original stage mismatch {ref['key']}: {stage} != {saved['earliest_stage']}")
            # Expanded span includes +/-6 snapping outside a source block.
            near = traces[bisect_left(ends, ref["span"][0] - 6):bisect_right(starts, ref["span"][1] + 6)]
            overlaps = [t for t in near if D.interval_overlap(ref["cds"], t["runs"]) > 0]
            candidates = [t for t in near if contains_truth(ref, t)]
            subruns = [t for t in near if proper_subruns(ref, t, sequence)]
            exact_chosen = chosen_map.get(chain(ref["cds"]), [])
            emitted = any(t.get("emitted_model") for t in exact_chosen)
            if emitted != (stage == "emitted_exact_chain"):
                raise AssertionError("exact emitted chains do not reconcile")
            category = primary_category(ref, emitted, exact_chosen, candidates, subruns, support, len(overlaps), cover[ref["key"]])
            exact_phase = [t for t in exact_chosen if t["failure_stage"] == "phase_check"]
            phase = phase_evidence(exact_phase[0], ref) if exact_phase else None
            failure_conflict = token_phase_conflict(ref, phase["failure"]["position"], sequence) if phase else None
            ref_start_conflicts = sum(token_phase_conflict(ref, a, sequence) for a, _b in ref["cds"])
            row = {"species": key[0], "seqid": key[1], "strand": key[2], "transcript_id": ref["transcript_id"],
                   "canonical": ref["canonical"], "CDS_count": len(ref["cds"]), "span_gt_6144": ref["span_gt_6144"],
                   "original_stage": stage, "category": category, "original_lineage_id": match["lineage_id"] if match else None,
                   "raw_covering_lineage_id": cover[ref["key"]]["lineage_id"] if cover[ref["key"]] else None,
                   "CDS_state_overlap_bp": support, "CDS_supporting_blocks": len(overlaps),
                   "truth_candidate_lineages": [t["lineage_id"] for t in candidates],
                   "proper_subrun_lineages": [t["lineage_id"] for t in subruns],
                   "exact_chosen_lineages": [t["lineage_id"] for t in exact_chosen],
                   "motif_reachable": int(saved["motif_reachable"]),
                   "phase_evidence": phase, "phase_failure_token_label_conflict": failure_conflict,
                   "truth_CDS_starts_with_token_label_conflict": ref_start_conflicts,
                   "annotation_phase_matches_cumulative": ref["annotation_phase_matches_cumulative"]}
            rows.append(row)
        del states
        print(f"epoch={epoch} reconciled {key}", flush=True)
    if len(rows) != 6450 or len({(r['species'], r['seqid'], r['strand'], r['transcript_id']) for r in rows}) != 6450:
        raise AssertionError("reference accounting failed")
    counts = Counter(r["original_stage"] for r in rows)
    if {s: counts[s] for s in D.STAGES} != old["reference_attrition"]["stage_counts"]:
        raise AssertionError("main stage totals changed")
    if len(phase_rows) != old["candidate_lineage_reconciliation"]["terminal_status_counts"]["phase_check"]:
        raise AssertionError("phase lineage terminal accounting failed")
    strata = defaultdict(Counter)
    for r in rows:
        for label in ("all", r["species"], "canonical" if r["canonical"] else "noncanonical",
                      "single_CDS" if r["CDS_count"] == 1 else "multi_CDS",
                      "span_gt_6144" if r["span_gt_6144"] else "span_lte_6144"):
            strata[label][r["category"]] += 1
    summary = {"epoch": epoch, "reference_total": len(rows), "exact": counts["emitted_exact_chain"],
               "missing": len(rows) - counts["emitted_exact_chain"], "category_counts_by_stratum": dict(strata),
               "original_stage_counts": dict(counts),
               "phase_lineage_total": len(phase_rows), "phase_lineage_evidence": dict(Counter(p["label"] for p in phase_rows)),
               "exact_phase_first_failure_CDS_index": dict(Counter(p["failed_CDS_index"] for p in phase_rows if p["reference_relation"] == "exact_chosen_chain")),
               "exact_phase_first_failure_expected_predicted": dict(Counter(f"{p['failure']['expected']}->{p['failure']['predicted']}" for p in phase_rows if p["reference_relation"] == "exact_chosen_chain")),
               "truth_start_token_conflicts": sum(r["truth_CDS_starts_with_token_label_conflict"] for r in rows),
               "truth_CDS_starts": sum(r["CDS_count"] for r in rows),
               "exact_phase_failure_token_conflicts": sum(r["phase_failure_token_label_conflict"] is True for r in rows),
               "exact_phase_failure_token_conflict_denominator": sum(r["phase_evidence"] is not None for r in rows),
               "annotation_phase_convention_disagreements": sum(not r["annotation_phase_matches_cumulative"] for r in rows),
               "recovered_truth_start_token_conflicts": sum(r["truth_CDS_starts_with_token_label_conflict"] for r in rows if r["category"] == "recovered_exact"),
               "recovered_truth_CDS_starts": sum(r["CDS_count"] for r in rows if r["category"] == "recovered_exact"),
               "phase_censoring": "Only the checked prefix and first failure are saved; downstream logits and score-filter outcomes are unidentifiable."}
    dst = output / f"epoch_{epoch}"
    dst.mkdir()
    for filename, records in (("reference_decomposition.jsonl", rows), ("phase_lineage_evidence.jsonl", phase_rows)):
        with (dst / filename).open("w") as h:
            for row in records:
                h.write(json.dumps(row, sort_keys=True) + "\n")
    D.save_json_atomic(dst / "summary.json", summary)
    return summary


def main():
    output = ROOT / "outputs" / OUT_ID
    if any(output.glob("epoch_*")) or (output / "summary.json").exists():
        raise FileExistsError("retrospective output already contains results")
    output.mkdir(exist_ok=True)
    sources = [ROOT / "outputs" / case.SOURCE_ID, ROOT / "outputs" / case.OUTPUT_ID]
    before = [case.source_stats(p) for p in sources]
    config = yaml.safe_load((ROOT / "configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml").read_text())
    species = D.load_species(ROOT, config)
    validation, _lengths = D.validation_truth(species)
    if set(validation) != case.VALIDATION:
        raise AssertionError("development scope mismatch")
    references = D.build_reference_records(species, validation)
    by_key = {(sp, sid, tx["strand"], tx["id"]): tx for (sp, sid), txs in validation.items() for tx in txs}
    for r in references:
        tx = D.m25._oriented_transcript(by_key[r["key"]], len(species[r["species"]]["seqs"][r["seqid"]]), r["strand"])
        codes = [int(p) + 1 if p in {"0", "1", "2"} else None for _a, _b, p in tx["CDS"]]
        r["phase_codes"] = codes
        r["phase_by_start"] = {a: p for (a, _b), p in zip(r["cds"], codes)}
        cumulative = 0
        consistent = True
        for (a, b), p in zip(r["cds"], codes):
            consistent &= p == 1 + (3 - cumulative % 3) % 3
            cumulative += b - a
        r["annotation_phase_matches_cumulative"] = consistent
    epochs = [analyze_epoch(e, species, references, output) for e in (1, 2, 3)]
    if before != [case.source_stats(p) for p in sources]:
        raise AssertionError("original R4 or correction files changed")
    D.save_json_atomic(output / "summary.json", {"status": "COMPLETED_RETROSPECTIVE_REVIEW_REQUIRED", "epochs": epochs,
                       "source_metadata_unchanged": True, "setaria_files_read": False, "model_inference": False,
                       "truth_assisted_total": "NOT_RECOVERABLE_WITH_SAVED_ARTIFACTS"})
    (output / "STATUS").write_text("COMPLETED_RETROSPECTIVE_REVIEW_REQUIRED\n")
    print("COMPLETED_RETROSPECTIVE_REVIEW_REQUIRED", flush=True)


if __name__ == "__main__":
    main()
