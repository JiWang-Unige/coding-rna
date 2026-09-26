import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/phase_panel_preflight.py"
spec = importlib.util.spec_from_file_location("phase_panel", PATH)
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


def test_proportional_rounding_is_deterministic_and_conserves_target():
    counts = {("a",): 7, ("b",): 2, ("c",): 1}
    assert P.allocate(counts, 6) == {("a",): 4, ("b",): 1, ("c",): 1}
    with pytest.raises(ValueError):
        P.allocate(counts, 11)


def test_selection_is_numeric_id_not_input_order_or_window_size():
    rows = [{"species": "a", "seqid": "c", "strand": "+", "lineage_id": i,
             "panel_group": "emitted_control", "block_span": [i * 1000000, i * 1000000 + 10]}
            for i in [10, 2, 1]]
    chosen, _ = P.choose(rows, 2)
    assert [r["lineage_id"] for r in chosen] == [1, 2]


def test_window_coverage_includes_competitors_and_acceptor_plus_two():
    trace = {"runs": [[5, 20], [6144, 6170]],
             "candidate_positions": {"start": [[5, 12290]], "stop": [[6167]], "donor": [[20]], "acceptor": [[6142]]},
             "chosen_positions": {"start": [5], "stop": [6167], "donor": [20], "acceptor": [6142]}}
    row = {"species": "a", "seqid": "c", "strand": "+", "lineage_id": 1, "trace": trace}
    positions = P.required_positions(trace)
    assert 6144 in positions and 12290 in positions
    assert [w["start"] for w in P.coverage([row])] == [0, 6144, 12288]


def test_orientations_do_not_share_inference_windows():
    trace = {"runs": [[5, 20]], "candidate_positions": {"start": [[5]], "stop": [[17]], "donor": [], "acceptor": []},
             "chosen_positions": {"start": [5], "stop": [17], "donor": [], "acceptor": []}}
    rows = [{"species": "a", "seqid": "c", "strand": strand, "lineage_id": 1, "trace": trace}
            for strand in ["+", "-"]]
    assert len(P.coverage(rows)) == 2


def test_window_budget_never_truncates_panel_or_coverage():
    rows = []
    for i in range(65):
        start = i * P.WINDOW + 5
        trace = {"runs": [[start, start + 15]],
                 "candidate_positions": {"start": [[start]], "stop": [[start + 12]], "donor": [], "acceptor": []},
                 "chosen_positions": {"start": [start], "stop": [start + 12], "donor": [], "acceptor": []}}
        rows.append({"species": "a", "seqid": "c", "strand": "+", "lineage_id": i, "trace": trace})
    assert len(P.coverage(rows)) == 65
    assert len(rows) == 65
