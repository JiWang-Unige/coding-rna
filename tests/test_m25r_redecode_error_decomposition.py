import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/redecode_error_decomposition.py"
SPEC = importlib.util.spec_from_file_location("m25r_redecode", SCRIPT)
DIAG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIAG)
M25 = DIAG.m25


def logits_for_states(states):
    logits = np.full((len(states), 3), -8.0, dtype=np.float32)
    logits[np.arange(len(states)), states] = 8.0
    return logits


def coding_example():
    sequence = list("A" * 40)
    sequence[5:12] = list("ATGAAAA")
    sequence[12:14] = list("GT")
    sequence[18:20] = list("AG")
    sequence[20:28] = list("AAAAATAA")
    states = np.array([M25.I] * 5 + [M25.C] * 7 + [M25.G] * 8 + [M25.C] * 8 + [M25.I] * 12)
    boundary_logits = np.full((40, 4), -20.0, dtype=np.float32)
    boundary_logits[5, 0] = 10.0
    boundary_logits[25, 1] = 10.0
    boundary_logits[12, 2] = 10.0
    boundary_logits[18, 3] = 10.0
    phase_logits = np.zeros((40, 4), dtype=np.float32)
    phase_logits[5, 1] = 10.0
    phase_logits[20, 3] = 10.0
    return "".join(sequence), logits_for_states(states), boundary_logits, phase_logits


def test_selects_one_original_highest_ranked_row_per_epoch():
    rows = []
    for epoch in (1, 2, 3):
        rows.extend([
            {"epoch": epoch, "rank": [0.1, 0.2, -0.01, -2], "enumeration_order": 2},
            {"epoch": epoch, "rank": [0.1, 0.2, -0.01, -1], "enumeration_order": 1},
            {"epoch": epoch, "rank": [0.2, 0.1, -0.02, -3], "enumeration_order": 3},
        ])
    selected = DIAG.select_epoch_rows({"rows": rows})
    assert {epoch: row["enumeration_order"] for epoch, row in selected.items()} == {1: 3, 2: 3, 3: 3}


def test_trace_reproduces_prefilter_and_frozen_filter_outputs():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    _states, traces, prefilter, emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    assert len(prefilter) == 1
    assert len(emitted) == 1
    assert traces[0]["failure_stage"] is None
    assert emitted[0]["cds"] == [(5, 12), (20, 28)]

    strict = dict(thresholds, donor=0.99999)
    _states, traces, prefilter, emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, strict
    )
    assert len(prefilter) == 1
    assert emitted == []
    assert traces[0]["failure_stage"] == "boundary_threshold_filter"


def test_reference_assignment_is_complete_and_uses_earliest_exact_failure():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    states, traces, _prefilter, _emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    reference = {
        "key": ("species", "chr", "+", "tx"),
        "species": "species",
        "seqid": "chr",
        "strand": "+",
        "transcript_id": "tx",
        "cds": [(5, 12), (20, 28)],
        "events": {"start": [5], "stop": [25], "donor": [12], "acceptor": [18]},
        "span": (5, 28),
        "canonical": True,
        "span_gt_6144": False,
        "tile_edge_within_6": True,
        "CDS_count": 2,
    }
    assignments, upper, _multiplicity, _matched = DIAG.assign_references(
        [reference], {("species", "chr", "+"): states},
        {("species", "chr", "+"): traces}, {("species", "chr", "+"): sequence},
        {("species", "chr", "+"): (boundary, phase)}, thresholds,
    )
    assert assignments[reference["key"]] == "emitted_exact_chain"
    assert upper == {
        "transition_reachable": 1,
        "motif_reachable": 1,
        "truth_assisted_exact_chain": 1,
    }

    changed = dict(reference, key=("species", "chr", "+", "shifted"), transcript_id="shifted")
    changed["events"] = dict(reference["events"], donor=[11])
    assignments, _upper, _multiplicity, _matched = DIAG.assign_references(
        [changed], {("species", "chr", "+"): states},
        {("species", "chr", "+"): traces}, {("species", "chr", "+"): sequence},
        {("species", "chr", "+"): (boundary, phase)}, thresholds,
    )
    assert assignments[changed["key"]] == "ordered_donor_acceptor_candidates"


def test_event_matching_is_one_to_one():
    assert DIAG.event_match_count([10, 12], [11], 1) == 1
    assert DIAG.event_match_count([10, 20], [9, 21], 1) == 2
    assert DIAG.event_match_count([10], [12], 1) == 0


def test_float16_histogram_average_precision_preserves_ties():
    scores = np.array([2.0, 1.0, 1.0, -1.0], dtype=np.float16)
    labels = np.array([1, 0, 1, 0], dtype=np.int64)
    keys = DIAG.float16_order_key(scores).astype(np.int64)
    total = np.bincount(keys, minlength=65536)
    positive = np.bincount(keys, weights=labels, minlength=65536).astype(np.int64)
    assert abs(DIAG.average_precision_from_hist(total, positive) - (0.5 + (2 / 3) * 0.5)) < 1e-12


def test_independent_gff_audit_checks_every_emitted_transcript():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    _states, _traces, _prefilter, emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    emitted[0]["strand"] = "+"
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "prediction.gff3"
        M25._write_gff3(str(path), {"chr": emitted}, "test")
        audit = DIAG.audit_gff3(path, {"chr": sequence}, {"chr": len(sequence)}, 1)
    assert audit["n_checked"] == 1
    assert audit["audit_coverage"] == 1.0
    assert audit["validity_fraction"] == 1.0
    assert audit["invalid_transcripts"] == 0
    assert audit["failures"] == []


def test_invalid_structure_is_recorded_without_aborting_complete_audit():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    _states, _traces, _prefilter, emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    emitted[0]["strand"] = "+"
    emitted[0]["phase"][1] = 0
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "prediction.gff3"
        M25._write_gff3(str(path), {"chr": emitted}, "test")
        audit = DIAG.audit_gff3(path, {"chr": sequence}, {"chr": len(sequence)}, 1)

    assert audit["audit_coverage"] == 1.0
    assert audit["valid_transcripts"] == 0
    assert audit["invalid_transcripts"] == 1
    assert audit["component_failure_counts"]["phase_continuity"] == 1
    assert audit["transcript_ledger"][0]["failed_components"] == ["phase_continuity"]


def test_gff_count_and_parent_linkage_mismatches_remain_fatal():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    _states, _traces, _prefilter, emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    emitted[0]["strand"] = "+"
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "prediction.gff3"
        M25._write_gff3(str(path), {"chr": emitted}, "test")
        with pytest.raises(AssertionError, match="model count"):
            DIAG.audit_gff3(path, {"chr": sequence}, {"chr": len(sequence)}, 2)

        text = path.read_text(encoding="utf-8")
        path.write_text(
            text.replace("Parent=M25_gene_0000001\n", "Parent=missing_gene\n", 1),
            encoding="utf-8",
        )
        with pytest.raises(AssertionError, match="linkage audit"):
            DIAG.audit_gff3(path, {"chr": sequence}, {"chr": len(sequence)}, 1)


def test_empty_gff_audit_is_complete_but_not_applicable():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "prediction.gff3"
        path.write_text("##gff-version 3\n", encoding="utf-8")
        audit = DIAG.audit_gff3(path, {"chr": "A"}, {"chr": 1}, 0)

    assert audit["complete_empty_audit"] is True
    assert audit["audit_coverage"] == 1.0
    assert audit["validity_fraction"] == "not_applicable"
    assert audit["transcript_ledger"] == []


def test_transition_upper_bound_does_not_reuse_one_anchor(monkeypatch):
    monkeypatch.setattr(
        DIAG.m25,
        "transition_candidates",
        lambda _states: [(11, M25.C, M25.G, ("donor",))],
    )
    sequence = list("A" * 30)
    sequence[10:12] = list("GT")
    sequence[12:14] = list("GT")
    reference = {"events": {"donor": [10, 12]}}
    assert DIAG.transition_upper_bounds(np.zeros(30), "".join(sequence), reference) == (False, False)


def test_motif_upper_bound_requires_each_truth_coordinate_in_candidate_union(monkeypatch):
    monkeypatch.setattr(
        DIAG.m25,
        "transition_candidates",
        lambda _states: [
            (10, M25.C, M25.G, ("donor",)),
            (12, M25.C, M25.G, ("donor",)),
        ],
    )
    sequence = list("A" * 30)
    sequence[12:14] = list("GT")
    reference = {"events": {"donor": [10, 12]}}
    assert DIAG.transition_upper_bounds(np.zeros(30), "".join(sequence), reference) == (True, False)


def test_truth_assisted_chain_still_obeys_phase_head():
    sequence, region, boundary, phase = coding_example()
    thresholds = {"region": 0.5, "start": 0.5, "stop": 0.5, "donor": 0.5, "acceptor": 0.5}
    _states, traces, _prefilter, _emitted = DIAG.trace_decode_orientation(
        sequence, region, boundary, phase, thresholds
    )
    reference = {
        "cds": [(5, 12), (20, 28)],
        "events": {"start": [5], "stop": [25], "donor": [12], "acceptor": [18]},
    }
    assert DIAG.truth_assisted_exact_chain(
        reference, traces[0], sequence, boundary, phase, thresholds
    )
    wrong_phase = phase.copy()
    wrong_phase[5] = 0.0
    wrong_phase[5, 2] = 10.0
    assert not DIAG.truth_assisted_exact_chain(
        reference, traces[0], sequence, boundary, wrong_phase, thresholds
    )


def test_main_persists_each_epoch_before_starting_the_next(monkeypatch):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "source"
        output = root / "output"
        source.mkdir()
        config_path = root / "config.yaml"
        config_path.write_text(json.dumps({
            "exp_id": "source",
            "seed": 0,
            "model": {"window_bp": DIAG.WINDOW_BP},
            "training": {
                "epochs": 3,
                "sample_fraction": 1.0,
                "train_window_cap": 0,
            },
            "data": {"development_species": [], "train_split": "train"},
        }), encoding="utf-8")
        (source / "validation_grid_diagnostics.json").write_text(json.dumps({
            "rows": [
                {"epoch": epoch, "rank": [1.0], "enumeration_order": epoch}
                for epoch in (1, 2, 3)
            ]
        }), encoding="utf-8")

        monkeypatch.setattr(DIAG, "EXPECTED_REFERENCE_CHAINS", 0)
        monkeypatch.setattr(DIAG, "EXPECTED_TRAIN_WINDOWS", 0)
        monkeypatch.setattr(DIAG, "load_species", lambda _root, _config: {})
        monkeypatch.setattr(DIAG, "validation_truth", lambda _species: ({}, {}))
        monkeypatch.setattr(DIAG, "build_reference_records", lambda _species, _refs: [])
        monkeypatch.setattr(
            DIAG,
            "load_inference_model",
            lambda _config: (None, None, None, None, None, 1),
        )
        monkeypatch.setattr(DIAG.m25, "OrientationWindowDataset", lambda *_args: [])

        def fake_run_epoch(epoch, *_args):
            out_dir = _args[-1]
            if epoch > 1:
                assert (out_dir / f"epoch_{epoch - 1}" / "diagnostic.json").is_file()
            result = {
                "reproduction": {"metric": {"absolute_error": 0.0}},
                "reference_attrition": {"stage_counts": {}},
                "independent_GFF3_validity": {
                    "audit_coverage": 1.0,
                    "invalid_transcripts": int(epoch == 1),
                },
            }
            DIAG.save_json_atomic(out_dir / f"epoch_{epoch}" / "diagnostic.json", result)
            return result

        monkeypatch.setattr(DIAG, "run_epoch", fake_run_epoch)
        monkeypatch.setattr(sys, "argv", [
            str(SCRIPT),
            "--root", str(root),
            "--experiment-id", "M25R-R3-TEST",
            "--config", str(config_path),
            "--m25r-output", str(source),
            "--out-dir", str(output),
        ])
        DIAG.main()

        assert all((output / f"epoch_{epoch}" / "diagnostic.json").is_file()
                   for epoch in (1, 2, 3))
        final = json.loads((output / "stage1_diagnostic.json").read_text(encoding="utf-8"))
        assert final["status"] == "COMPLETED_STAGE1_REVIEW_REQUIRED"
        assert final["scientific_status"] == "SCIENTIFIC_NO_GO_INVALID_STRUCTURES"
        assert (output / "STATUS").read_text(encoding="utf-8").strip() == final["status"]
