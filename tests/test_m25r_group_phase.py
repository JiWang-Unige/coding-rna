import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("group_phase", ROOT / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/group_phase_retrospective.py")
G = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(G)


def reference(cds=((5, 12), (20, 28))):
    return {"key": ("sp", "chr", "+", "tx"), "transcript_id": "tx", "cds": list(cds),
            "events": G.D.truth_events(cds), "span": (cds[0][0], cds[-1][1]),
            "phase_by_start": {5: 1, 20: 3}, "phase_codes": [1, 3]}


def trace(cds=((5, 12), (20, 28)), predicted=1):
    return {"lineage_id": 1, "block_span": [3, 30], "runs": list(cds),
            "chosen_positions": G.D.truth_events(cds),
            "candidate_positions": {k: [[p] for p in v] for k, v in G.D.truth_events(cds).items()},
            "phase_checks": [{"position": 5, "expected": 1, "predicted": 1},
                             {"position": 20, "expected": 3, "predicted": predicted}],
            "failure_stage": "phase_check"}


def test_correct_chain_phase_failure_is_not_grouping_failure():
    r, t = reference(), trace()
    assert G.contains_truth(r, t)
    assert G.original_stage(r, t, 15) == "phase_check"
    assert G.phase_evidence(t, r)["label"] == "exact_chain_correct_start_phase_disagreement"
    assert G.primary_category(r, False, [t], [t], [], 15, 1, t) == "correct_chosen_chain_phase_rejection"


def test_wrong_boundary_shifts_expected_frame_despite_correct_local_prediction():
    r = reference()
    t = trace(((5, 13), (20, 28)), predicted=3)
    t["phase_checks"][-1]["expected"] = 2
    assert G.original_stage(r, t, 15) == "ordered_donor_acceptor_candidates"
    assert G.phase_evidence(t, r)["label"] == "wrong_chain_frame_shift_head_agrees_reference"


def test_wrong_boundary_at_nonreference_start_is_not_phase_accuracy():
    t = trace(((5, 12), (21, 28)), predicted=1)
    t["phase_checks"][-1]["position"] = 21
    assert G.phase_evidence(t, reference())["label"] == "off_reference_CDS_start_unidentifiable_phase"


def test_phase_censoring_requires_prefix_then_first_failure():
    t = trace()
    t["phase_checks"][0]["predicted"] = 2
    with pytest.raises(AssertionError, match="first rejection"):
        G.validate_phase_trace(t)


def test_matching_reserves_exact_emitted_model_before_overlap_tie():
    r = reference()
    other = dict(r, key=("sp", "chr", "+", "a_other"), transcript_id="a_other", cds=[(4, 12), (20, 28)])
    t = trace()
    t["failure_stage"] = None
    t["emitted_model"] = {"cds": r["cds"]}
    matched, _cover = G.original_matching([other, r], [t])
    assert matched == {r["key"]: t}


def test_compatible_subruns_are_candidate_evidence_not_phase_rescue():
    seq = list("A" * 50)
    seq[5:8] = "ATG"
    seq[12:14] = "GT"
    seq[18:20] = "AG"
    seq[25:28] = "TAA"
    t = trace()
    t["runs"] = [[5, 12], [20, 28], [40, 46]]
    assert G.proper_subruns(reference(), t, "".join(seq))
    assert G.primary_category(reference(), False, [], [], [t], 15, 1, t) == "whole_block_excludes_truth_compatible_subruns"


def test_phase_1x1_identical_inputs_can_have_conflicting_targets():
    r = {"cds": [(0, 6)], "phase_codes": [1]}
    assert G.token_phase_conflict(r, 0, "AAAAAA")
    assert not G.token_phase_conflict(r, 0, "ACGACG")


def test_annotation_phase_conflict_not_blindly_called_head_error():
    r = reference()
    r["phase_by_start"][20] = 2
    assert G.phase_evidence(trace(), r)["label"] == "exact_chain_annotation_decoder_phase_disagreement"
