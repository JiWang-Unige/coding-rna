#!/usr/bin/env python3
"""R2: one preapproved systematic panel; unchanged epoch-1 accounting, no GPU."""
from collections import defaultdict
import phase_panel_preflight as P

OUTPUT = P.ROOT / "outputs/M25R-E1-PHASE-GATE-PANEL-R2"


def midpoint_indices(total, selected):
    if not 0 <= selected <= total:
        raise ValueError("insufficient eligible rows for fixed quota")
    return [(2 * j + 1) * total // (2 * selected) for j in range(selected)]


def choose_systematic(pool, target):
    buckets = defaultdict(list)
    for row in pool:
        buckets[P.stratum(row)].append(row)
    quotas = P.allocate({key: len(rows) for key, rows in buckets.items()}, target)
    result, inventory = [], []
    for key in sorted(buckets):
        rows = sorted(buckets[key], key=P.identity)
        indices = midpoint_indices(len(rows), quotas[key])
        for index in indices:
            row = dict(rows[index])
            trace = row["trace"]
            positions = P.required_positions(trace)
            row["eligible_stratum_index_zero_based"] = index
            row["CDS_count"] = len(trace["runs"])
            row["saved_phase_check_count"] = len(trace["phase_checks"])
            row["chosen_chain_span_bp"] = trace["chosen_positions"]["stop"][0] + 3 - trace["chosen_positions"]["start"][0]
            row["score_position_within_6bp_of_tile_edge"] = any(
                min(p % P.WINDOW, P.WINDOW - 1 - p % P.WINDOW) <= 6 for p in positions)
            result.append(row)
        inventory.append({"stratum": key, "eligible": len(rows), "selected": quotas[key],
                          "indices_zero_based": indices})
    return result, inventory


if __name__ == "__main__":
    P.main(output=OUTPUT, selector=choose_systematic,
           selection="proportional_largest_remainder_then_stratum_midpoints_floor((j+0.5)*N/n)")
