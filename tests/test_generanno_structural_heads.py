import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
TRAINER_PATH = ROOT / "src/foundation_probe/train_generanno_structural_heads.py"
TRAINER_SPEC = importlib.util.spec_from_file_location("m25_trainer", TRAINER_PATH)
M25 = importlib.util.module_from_spec(TRAINER_SPEC)
TRAINER_SPEC.loader.exec_module(M25)

EVALUATOR_PATH = ROOT / "scripts/eval_m25_structure.py"
EVALUATOR_SPEC = importlib.util.spec_from_file_location("m25_evaluator", EVALUATOR_PATH)
EVAL = importlib.util.module_from_spec(EVALUATOR_SPEC)
EVALUATOR_SPEC.loader.exec_module(EVAL)


def transcript(strand, cds, exons=None, partial=False):
    return {
        "id": f"tx_{strand}",
        "strand": strand,
        "CDS": cds,
        "exon": exons or cds,
        "partial": partial,
    }


def logits_for_states(states):
    logits = np.full((len(states), 3), -8.0, dtype=np.float32)
    logits[np.arange(len(states)), states] = 8.0
    return logits


def test_strand_targets_exclude_opposite_orientation_and_freeze_coordinates():
    plus = transcript("+", [(10, 19, "0"), (30, 39, "0")])
    minus = transcript("-", [(10, 19, "0"), (30, 39, "0")])
    targets = M25.build_strand_targets(60, [plus, minus], "+")

    assert np.all(targets["region"][10:19] == M25.C)
    assert np.all(targets["region"][30:39] == M25.C)
    assert targets["boundary"][10, 0] == 1       # ATG starts at 0-based 10; GFF start 11.
    assert targets["boundary"][36, 1] == 1       # Stop codon starts at 36; GFF 37..39.
    assert targets["boundary"][19, 2] == 1       # Donor GT starts at 19; GFF 20..21.
    assert targets["boundary"][28, 3] == 1       # Acceptor AG starts at 28; GFF 29..30.
    assert targets["boundary"].sum() == 4

    minus_only = M25.build_strand_targets(60, [plus], "-")
    assert np.all(minus_only["region"] == M25.I)
    assert minus_only["boundary"].sum() == 0
    assert np.all(minus_only["phase"] == 0)


def test_reverse_complement_interval_mapping_is_involutive():
    assert M25.reverse_complement("ACGTN") == "NACGT"
    assert M25.map_rc_interval(10, 19, 60) == (41, 50)
    assert M25.map_rc_interval(*M25.map_rc_interval(10, 19, 60), 60) == (10, 19)

    model = {
        "cds": [(5, 12), (20, 28)],
        "phase": [0, 2],
        "start_codon": (5, 8),
        "stop_codon": (25, 28),
    }
    mapped = M25._map_model_from_orientation(model, 60, "-")
    assert mapped["strand"] == "-"
    assert mapped["cds"] == [(32, 40), (48, 55)]
    assert mapped["phase"] == [2, 0]
    assert mapped["start_codon"] == (52, 55)
    assert mapped["stop_codon"] == (32, 35)


def test_primary_chromosome_allowlist_is_exact():
    sequences = {"chr1": "AAAA", "chloroplast": "CCCC", "scaffold": "GGGG"}
    assert M25.select_primary_chromosomes(sequences, ["chr1"]) == {"chr1": "AAAA"}
    try:
        M25.select_primary_chromosomes(sequences, ["chr2"])
    except ValueError as error:
        assert "absent from FASTA" in str(error)
    else:
        raise AssertionError("missing primary chromosome was accepted")


def test_frozen_region_states_and_transition_candidates():
    states = np.array([M25.I, M25.G, M25.C, M25.G, M25.I, M25.C, M25.I])
    transitions = M25.transition_candidates(states)
    assert transitions == [
        (2, M25.G, M25.C, ("start", "acceptor")),
        (3, M25.C, M25.G, ("stop", "donor")),
        (5, M25.I, M25.C, ("start",)),
        (6, M25.C, M25.I, ("stop",)),
    ]

    decoded = M25.region_state_path(logits_for_states(states), 0.5)
    assert decoded.tolist() == states.tolist()


def test_phase_continuity_and_unchanged_input_ablation():
    sequence = list("A" * 40)
    sequence[5:12] = list("ATGAAAA")
    sequence[12:14] = list("GT")
    sequence[18:20] = list("AG")
    sequence[20:28] = list("AAAAATAA")
    sequence = "".join(sequence)
    states = np.array([M25.I] * 5 + [M25.C] * 7 + [M25.G] * 8 + [M25.C] * 8 + [M25.I] * 12)
    region_logits = logits_for_states(states)
    boundary_logits = np.full((40, 4), -20.0, dtype=np.float32)
    phase_logits = np.zeros((40, 4), dtype=np.float32)
    phase_logits[5, 1] = 10.0
    phase_logits[20, 3] = 10.0
    thresholds = {"region": 0.5, "start": 0.99, "stop": 0.99, "donor": 0.99, "acceptor": 0.99}

    assert M25.decode_orientation(
        sequence, region_logits, boundary_logits, phase_logits, thresholds
    ) == []
    ablation = M25.decode_orientation(
        sequence, region_logits, boundary_logits, phase_logits, thresholds, ablation=True
    )
    assert [(row["cds"], row["phase"]) for row in ablation] == [
        ([(5, 12), (20, 28)], [0, 2])
    ]

    boundary_logits[5, 0] = 10.0
    boundary_logits[25, 1] = 10.0
    boundary_logits[12, 2] = 10.0
    boundary_logits[18, 3] = 10.0
    full = M25.decode_orientation(
        sequence, region_logits, boundary_logits, phase_logits, thresholds
    )
    assert [(row["cds"], row["phase"]) for row in full] == [
        ([(5, 12), (20, 28)], [0, 2])
    ]

    phase_logits[20] = 0.0
    phase_logits[20, 2] = 10.0
    assert M25.decode_orientation(
        sequence, region_logits, boundary_logits, phase_logits, thresholds
    ) == []


def test_embargo_marker_hashes_are_verified_before_reference_use(tmp_path):
    paths = {}
    for key in EVAL.MARKER_KEYS:
        path = tmp_path / key
        path.write_text(key)
        paths[key] = path
    marker = tmp_path / "SETARIA_EMBARGO_RELEASED.json"
    marker.write_text(json.dumps({
        key: {"path": str(path), "sha256": EVAL.sha256(path)}
        for key, path in paths.items()
    }))
    cli = {
        "allowlist": paths["primary_chromosome_allowlist"].resolve(),
        "genome_fasta": paths["genome_fasta"].resolve(),
        "full_gff3": paths["full_prediction_gff3"].resolve(),
        "ablation_gff3": paths["ablation_prediction_gff3"].resolve(),
    }
    assert set(EVAL.verify_embargo(marker, cli)) == set(EVAL.MARKER_KEYS)

    paths["checkpoint"].write_text("changed")
    try:
        EVAL.verify_embargo(marker, cli)
    except ValueError as error:
        assert "sha256 mismatch" in str(error)
    else:
        raise AssertionError("modified frozen checkpoint was accepted")
