import importlib.util
import sys
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION"
sys.path.insert(0, str(SOURCE))
spec = importlib.util.spec_from_file_location("systematic_panel", SOURCE / "phase_panel_systematic.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def test_midpoints_cover_ordered_list_and_handle_zero_full_and_insufficient():
    assert S.midpoint_indices(100, 4) == [12, 37, 62, 87]
    assert S.midpoint_indices(7, 3) == [1, 3, 5]
    assert S.midpoint_indices(5, 5) == [0, 1, 2, 3, 4]
    assert S.midpoint_indices(5, 0) == []
    with pytest.raises(ValueError):
        S.midpoint_indices(2, 3)


def test_same_quota_unique_midpoints_not_prefix_or_score_selection():
    rows = []
    for i in reversed(range(10)):
        start = i * 6144 + 5
        trace = {"runs": [[start, start + 15]], "phase_checks": [{"position": start}],
                 "candidate_positions": {"start": [[start]], "stop": [[start + 12]], "donor": [], "acceptor": []},
                 "chosen_positions": {"start": [start], "stop": [start + 12], "donor": [], "acceptor": []}}
        rows.append({"species": "a", "seqid": "c", "strand": "+", "lineage_id": i,
                     "panel_group": "emitted_control", "trace": trace})
    chosen, inventory = S.choose_systematic(rows, 4)
    assert [r["lineage_id"] for r in chosen] == [1, 3, 6, 8]
    assert inventory[0]["selected"] == 4
    assert all(r["CDS_count"] == 1 and r["saved_phase_check_count"] == 1 for r in chosen)
    assert all(r["score_position_within_6bp_of_tile_edge"] for r in chosen)
    assert all("eligible_stratum_index_zero_based" not in r for r in rows)


def test_all_midpoint_indices_are_unique_and_in_range():
    for total in (1, 7, 49, 191, 716, 801):
        for quota in range(min(total, 64) + 1):
            indices = S.midpoint_indices(total, quota)
            assert len(indices) == len(set(indices)) == quota
            assert all(0 <= i < total for i in indices)
