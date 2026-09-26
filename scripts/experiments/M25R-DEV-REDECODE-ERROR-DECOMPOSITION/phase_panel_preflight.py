#!/usr/bin/env python3
"""Freeze an epoch-1 lineage panel before score access; count original windows.

Standard-library only. No FASTA, checkpoint, model, or Setaria access.
Selection is the previously proposed stratified numeric-lineage-ID prefix;
it is a diagnostic convenience panel, not a random population sample.
"""
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4"
RETRO = ROOT / "outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE"
OUT = ROOT / "outputs/M25R-E1-PHASE-GATE-PANEL"
WINDOW = 6144
LIMIT = 64
GROUPS = {"exact_phase_failure": 64, "nonexact_phase_failure": 32, "emitted_control": 32}
SCOPE = {("arabidopsis_thaliana", "NC_003074.8"), ("oryza_sativa", "NC_089041.1")}


def read_jsonl(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle]


def identity(row):
    return row["species"], row["seqid"], row["strand"], row["lineage_id"]


def stratum(row):
    key = (row["species"], row["seqid"], row["strand"])
    if row["panel_group"] == "exact_phase_failure":
        key += ("conflict" if row["phase_failure_token_label_conflict"] else "no_conflict",)
    return key


def allocate(counts, target):
    """Proportional largest remainder; stable tuple ties, no score/size input."""
    total = sum(counts.values())
    if total < target:
        raise ValueError("insufficient eligible lineages; no substitute sampling")
    quotas = {key: target * count // total for key, count in counts.items()}
    left = target - sum(quotas.values())
    order = sorted(counts, key=lambda key: (-(target * counts[key] % total), key))
    for key in order[:left]:
        quotas[key] += 1
    return quotas


def choose(pool, target):
    buckets = defaultdict(list)
    for row in pool:
        buckets[stratum(row)].append(row)
    quotas = allocate({key: len(rows) for key, rows in buckets.items()}, target)
    result = []
    inventory = []
    for key in sorted(buckets):
        rows = sorted(buckets[key], key=lambda row: row["lineage_id"])
        result.extend(rows[:quotas[key]])
        inventory.append({"stratum": key, "eligible": len(rows), "selected": quotas[key]})
    return result, inventory


def required_positions(trace):
    """All motif-choice competitors plus every chosen CDS start for phase.

    Frozen downstream boundary scores are minima of chosen motif scores only.
    ORF sequence checks need FASTA but no further GPU positions. Region runs
    are held fixed: this is a candidate-conditional, not whole-pipeline replay.
    """
    positions = set()
    for rows in trace["candidate_positions"].values():
        for options in rows:
            positions.update(options)
    chosen = trace["chosen_positions"]
    for values in chosen.values():
        positions.update(values)
    if len(chosen["start"]) != 1 or len(chosen["stop"]) != 1:
        raise AssertionError("panel trace lacks a complete chosen chain")
    if len(chosen["acceptor"]) != len(trace["runs"]) - 1 or len(chosen["donor"]) != len(trace["runs"]) - 1:
        raise AssertionError("panel trace has incomplete splice choices")
    positions.update(p + 2 for p in chosen["acceptor"])
    if not positions or min(positions) < 0:
        raise AssertionError("invalid required coordinate")
    return sorted(positions)


def coverage(panel):
    users = defaultdict(list)
    for row in panel:
        positions = required_positions(row["trace"])
        windows = sorted({p // WINDOW * WINDOW for p in positions})
        row["required_score_positions"] = positions
        row["original_window_starts"] = windows
        for start in windows:
            users[(row["species"], row["seqid"], row["strand"], start)].append(identity(row))
    return [{"species": key[0], "seqid": key[1], "strand": key[2], "start": key[3],
             "window_bp": WINDOW, "lineages": users[key]} for key in sorted(users)]


def write_json(path, content):
    with path.open("x") as handle:
        json.dump(content, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main(output=OUT, selector=choose,
         selection="proportional_largest_remainder_then_numeric_lineage_id_prefix"):
    # R2 reuses the same accounting, with an explicitly separate selector/output.
    OUT = output
    if (OUT / "frozen_panel.json").exists() or (OUT / "preflight.json").exists():
        raise FileExistsError("panel or preflight already exists; no reselection")
    references = read_jsonl(RETRO / "epoch_1/reference_decomposition.jsonl")
    evidence = {identity(row): row for row in read_jsonl(RETRO / "epoch_1/phase_lineage_evidence.jsonl")}
    traces = read_jsonl(SOURCE / "epoch_1/candidate_lineages.jsonl")
    exact = {}
    for row in references:
        for lid in row["exact_chosen_lineages"]:
            key = row["species"], row["seqid"], row["strand"], lid
            if key in exact:
                raise AssertionError("ambiguous exact-chain reference")
            exact[key] = row
    pools = defaultdict(list)
    for trace in traces:
        key = identity(trace)
        if key[:2] not in SCOPE or key[2] not in {"+", "-"}:
            raise AssertionError("trace outside approved development scope")
        ref = exact.get(key)
        if trace["failure_stage"] == "phase_check":
            group = "exact_phase_failure" if ref else "nonexact_phase_failure"
            if key not in evidence:
                raise AssertionError("missing phase-lineage evidence")
        elif trace.get("emitted_model"):
            group = "emitted_control"
        else:
            continue
        conflict = ref["phase_failure_token_label_conflict"] if ref else None
        if group == "exact_phase_failure" and not isinstance(conflict, bool):
            raise AssertionError("missing exact first-failure conflict label")
        pools[group].append({"species": key[0], "seqid": key[1], "strand": key[2], "lineage_id": key[3],
                             "panel_group": group, "exact_reference_id": ref["transcript_id"] if ref else None,
                             "is_exact_to_primary_reference": ref is not None,
                             "phase_failure_token_label_conflict": conflict,
                             "phase_evidence": evidence.get(key), "trace": trace})
    expected = {"exact_phase_failure": 1149, "nonexact_phase_failure": 2129, "emitted_control": 2098}
    if {name: len(pools[name]) for name in GROUPS} != expected:
        raise AssertionError("eligible pool does not reconcile with frozen epoch-1 result")
    panel, inventory = [], {}
    for group, target in GROUPS.items():
        chosen, strata = selector(pools[group], target)
        panel.extend(chosen)
        inventory[group] = strata
    if len(panel) != 128 or len({identity(row) for row in panel}) != 128:
        raise AssertionError("panel must contain exactly 128 unique lineages")
    # Persist selected identities BEFORE evaluating the window cap. No replacement.
    OUT.mkdir(exist_ok=True)
    with (SOURCE / "epoch_1/diagnostic.json").open() as handle:
        frozen_tuple = json.load(handle)["frozen_tuple"]
    if frozen_tuple["epoch"] != 1 or frozen_tuple["enumeration_order"] != 601:
        raise AssertionError("wrong frozen epoch/tuple")
    write_json(OUT / "frozen_panel.json", {"selection": selection,
               "selection_uses_scores_or_window_coverage": False, "inference_performed": False,
               "strata": inventory, "frozen_tuple": frozen_tuple, "panel": panel})
    windows = coverage(panel)
    write_json(OUT / "window_coverage.json", {"windows": windows, "panel": panel})
    by_group = {group: len({(r["species"], r["seqid"], r["strand"], s)
                            for r in panel if r["panel_group"] == group for s in r["original_window_starts"]})
                for group in GROUPS}
    status = "WINDOW_PREFLIGHT_PASS" if len(windows) <= LIMIT else "STOP_WINDOW_BUDGET_EXCEEDED"
    summary = {"status": status, "panel_size": len(panel), "groups": dict(Counter(r["panel_group"] for r in panel)),
               "eligible_counts": expected, "window_count": len(windows), "window_limit": LIMIT,
               "window_count_by_group_nonadditive": by_group,
               "window_count_by_orientation": dict(Counter(f"{w['species']}|{w['strand']}" for w in windows)),
               "gpu_job_submitted": False, "gpu_inference_performed": False,
               "setaria_files_read": False, "weights_updated": False,
               "population_inference": False, "reference_thresholds_are_exploratory_only": True,
               "original_results_modified": False}
    write_json(OUT / "preflight.json", summary)
    (OUT / "STATUS").write_text(status + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
