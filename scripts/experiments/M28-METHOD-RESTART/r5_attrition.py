#!/usr/bin/env python3
"""Read-only R5 attrition diagnosis; never changes candidates, scores or calls."""
import argparse
import bisect
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import common_ruler as E
from src.m28.assembly import nearest_owner
from src.m28.c0_decode import genomic_chain

def chain_key(tx):
    return tuple((a, b) for a, b, *_ in tx["CDS"])

def quantiles(values):
    values = sorted(values)
    if not values:
        return {"n": 0}
    return {"n": len(values), **{
        name: values[round((len(values)-1)*p)]
        for name, p in (("min", 0), ("p50", .5), ("p90", .9), ("p99", .99), ("max", 1))
    }}

def geometry(chain):
    return {"exons": len(chain), "span_bp": chain[-1][1]-chain[0][0],
            "CDS_bp": sum(b-a for a, b in chain),
            "max_intron_bp": max((v[0]-u[1] for u, v in zip(chain, chain[1:])), default=0)}

def main(pilot, out):
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    result = json.loads((pilot/"evaluation/summary.json").read_text())
    infer = json.loads((pilot/"infer_B1/summary.json").read_text())
    refs = {}
    lengths = {seqid: length for _, seqid, length, _ in E.S.SCOPE}
    for species, seqid, _, _ in E.S.SCOPE:
        with (E.ROOT/"outputs/M28-COMMON-RULER-R3/result"/(species+".references.jsonl")).open() as f:
            refs[species] = [json.loads(line) for line in f]
    predictions = {}
    for method in ("C0", "B1"):
        gff = pilot/"evaluation"/(method+".gff3")
        predictions[method] = E.S.E.primary_transcripts(
            E.S.E.parse_annotation(gff, lengths), complete_only=False)
    report = {"experiment": "M28-PILOT-R5-ATTRITION-D1",
              "status": "exploratory_read_only_frozen_outputs",
              "parent": "M28-PILOT-R5", "new_forward_pass": False,
              "new_training": False, "prediction_or_threshold_changed": False,
              "test_setaria_used": False, "scopes": {}, "per_species": {}}
    all_cases = []
    for scope in infer["scopes"]:
        species, seqid, strand = (scope[k] for k in ("species", "seqid", "strand"))
        name = species+"|"+strand
        windows = scope["windows"]
        starts, ends = [w[0] for w in windows], [w[1] for w in windows]
        centers = [a+b for a, b in windows]
        selected = {chain_key(t) for t in predictions["B1"] if t["seqid"] == seqid and t["strand"] == strand}
        c0 = {chain_key(t) for t in predictions["C0"] if t["seqid"] == seqid and t["strand"] == strand}
        reference = {chain_key(t): t for t in refs[species] if t["strand"] == strand}
        candidates = {c: [] for c in reference}
        owned = {}
        contained = defaultdict(list)
        support = {}
        for c in reference:
            indices = list(range(bisect.bisect_left(ends, c[-1][1]),
                                 bisect.bisect_right(starts, c[0][0])))
            owner = nearest_owner(c, windows, centers)
            support[c] = {"containing_windows": len(indices), "owner_index": owner,
                          "owner_contains": owner in indices,
                          "all_links_in_any_containing_window": False,
                          "all_links_in_owner": False}
            for i in indices:
                contained[i].append(c)
        counts = Counter()
        finite = True
        with (pilot/"infer_B1"/scope["file"]).open() as f:
            for i, line in enumerate(f):
                record = json.loads(line)
                wa, wb = windows[i]
                expected_id = species+"|"+seqid+"|"+str(wa)+"|"+strand
                if record["window_index"] != i or record["window_id"] != expected_id or record["reference_used"] is not False:
                    raise ValueError("Frozen window identity/provenance mismatch")
                if len(record["chains"]) != len(record["score"]):
                    raise ValueError("Candidate/score length mismatch")
                counts["windows"] += 1
                counts["window_candidate_records"] += len(record["chains"])
                budget = record["budget"]
                cap = max(1, math.ceil((wb-wa)/1024)*budget["chains_per_kb"])
                counts["windows_at_final_cap"] += len(record["chains"]) == cap
                for field in ("complete_paths", "final_chains_pruned", "beam_states_pruned",
                              "intron_links_considered", "intron_links_retained"):
                    counts[field] += record["counts"][field]
                links = {tuple(pair) for pair in record["links"]}
                for c in contained[i]:
                    oriented = ([(a-wa, b-wa) for a, b in c] if strand == "+" else
                                [(wb-b, wb-a) for a, b in reversed(c)])
                    required = {(u[1], v[0]) for u, v in zip(oriented, oriented[1:])}
                    if required <= links:
                        support[c]["all_links_in_any_containing_window"] = True
                        if support[c]["owner_index"] == i:
                            support[c]["all_links_in_owner"] = True
                for chain, score in zip(record["chains"], record["score"]):
                    if not math.isfinite(score):
                        finite = False
                    c = tuple(genomic_chain(chain, wa, wb, strand))
                    if c not in reference:
                        continue
                    candidates[c].append((i, float(score)))
                    if support[c]["owner_index"] == i:
                        owned[c] = max(float(score), owned.get(c, -math.inf))
        if counts["windows"] != len(windows) or not finite:
            raise ValueError("Partial inference or nonfinite score")
        observed = {"free_unique": sum(bool(v) for v in candidates.values()),
                    "owner_unique": len(owned), "selected_unique": len(selected & reference.keys())}
        if observed != result["methods"]["B1"]["stages"][name]["exact_reference_hits"]:
            raise ValueError(("Attrition counts disagree with frozen evaluation", name, observed))
        if counts["window_candidate_records"] != scope["candidate_records"]:
            raise ValueError("Candidate record count differs from inference summary")
        categories = Counter()
        for c, tx in reference.items():
            seen = candidates[c]
            if c in selected:
                category = "selected_exact"
            elif c in owned:
                category = "owned_nonpositive" if owned[c] <= 0 else "owned_positive_conflict"
            elif seen:
                category = "free_nonowner_only"
            else:
                category = "absent_from_final_free"
            if c in selected and (c not in owned or owned[c] <= 0):
                raise ValueError("Selected exact chain is not a positive owner candidate")
            overlaps = [p for p in selected if p[0][0] < c[-1][1] and p[-1][1] > c[0][0]]
            case = {"species": species, "seqid": seqid, "strand": strand,
                    "reference_id": tx["id"], "category": category, **geometry(c), **support[c],
                    "unsupported_labels": tx["M28_unsupported_strata_not_excluded"],
                    "C0_selected_exact": c in c0, "seen_window_count": len(seen),
                    "best_observed_gain": max((v for _, v in seen), default=None),
                    "owner_gain": owned.get(c),
                    "selected_overlap_count": len(overlaps),
                    "overlap_with_other_exact_reference": any(p in reference and p != c for p in overlaps),
                    "overlap_with_reference_nonexact": any(p not in reference for p in overlaps)}
            all_cases.append(case)
            categories[category] += 1
        report["scopes"][name] = {"reference": len(reference), "replayed_hits": observed,
                                 "category_counts": dict(categories), "proposal_counts": dict(counts)}
    if len(all_cases) != result["reference_primary_chains"]:
        raise ValueError("Primary reference denominator mismatch")
    for species in refs:
        cases = [c for c in all_cases if c["species"] == species]
        categories = Counter(c["category"] for c in cases)
        strata = {}
        for label, get in (
            ("exon_count", lambda c: "1" if c["exons"] == 1 else "2-4" if c["exons"] <= 4 else "5-8" if c["exons"] <= 8 else "9+"),
            ("span", lambda c: "<=1kb" if c["span_bp"] <= 1000 else "1-5kb" if c["span_bp"] <= 5000 else "5-12kb" if c["span_bp"] <= 12000 else "12-24.576kb" if c["span_bp"] <= 24576 else ">24.576kb")):
            strata[label] = {}
            for c in cases:
                row = strata[label].setdefault(get(c), Counter())
                row["reference"] += 1
                row[c["category"]] += 1
                row["C0_selected_exact"] += c["C0_selected_exact"]
        absent = [c for c in cases if c["category"] == "absent_from_final_free"]
        nonowner = [c for c in cases if c["category"] == "free_nonowner_only"]
        conflict = [c for c in cases if c["category"] == "owned_positive_conflict"]
        prediction_geometry = {}
        ref_keys = {(t["strand"], chain_key(t)) for t in refs[species]}
        for method in ("C0", "B1"):
            for exact in (True, False):
                chains = [chain_key(t) for t in predictions[method]
                          if t["seqid"] == next(q for s, q, _, _ in E.S.SCOPE if s == species)
                          and ((t["strand"], chain_key(t)) in ref_keys) == exact]
                vals = [geometry(c) for c in chains]
                prediction_geometry[method+("_exact" if exact else "_nonexact")] = {
                    k: quantiles([v[k] for v in vals]) for k in ("exons", "span_bp", "CDS_bp", "max_intron_bp")}
        report["per_species"][species] = {
            "reference": len(cases), "category_counts": dict(categories), "strata": strata,
            "absent_breakdown": {
                "not_contained_in_any_window": sum(c["containing_windows"] == 0 for c in absent),
                "contained_but_necessary_retained_links_absent_in_every_window": sum(c["containing_windows"] > 0 and not c["all_links_in_any_containing_window"] for c in absent),
                "retained_links_available_but_final_chain_absent": sum(c["all_links_in_any_containing_window"] for c in absent)},
            "nonowner_diagnostic": {
                "owner_does_not_contain": sum(not c["owner_contains"] for c in nonowner),
                "owner_necessary_links_absent": sum(not c["all_links_in_owner"] for c in nonowner),
                "positive_gain_in_some_nonowner_window": sum(c["best_observed_gain"] > 0 for c in nonowner)},
            "positive_conflict_diagnostic": {
                "with_reference_nonexact_selected": sum(c["overlap_with_reference_nonexact"] for c in conflict),
                "with_other_exact_selected": sum(c["overlap_with_other_exact_reference"] for c in conflict),
                "without_selected_overlap": sum(c["selected_overlap_count"] == 0 for c in conflict)},
            "prediction_geometry": prediction_geometry}
    report["limits"] = [
        "DEV-selected exploratory description, not independent confirmation or a new method score.",
        "Missing retained links are a necessary-support failure, not isolated link-network causation: endpoint pruning also removes links.",
        "Retained links do not establish exon legality, phase-compatible paths, or survival of exon/beam/final pruning.",
        "Nonowner gains are observed values, not calibrated scores transferable across windows.",
        "Overlap counts are nonexclusive and do not prove biological false positives.",
        "No decoder replay, no new rank table, no parameter or candidate-budget search."]
    with (out/"reference_cases.jsonl").open("x") as f:
        for row in all_cases:
            f.write(json.dumps(row)+"\n")
    with (out/"summary.json").open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps({"status": report["status"], "per_species": {
        s: {k: v for k, v in r.items() if k != "prediction_geometry"}
        for s, r in report["per_species"].items()}}), flush=True)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pilot-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    main(a.pilot_dir.resolve(), a.output_dir.resolve())
