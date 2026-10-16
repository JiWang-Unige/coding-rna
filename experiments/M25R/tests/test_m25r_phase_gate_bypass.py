import importlib.util
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

SOURCE = Path(__file__).resolve().parents[1] / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION"
sys.path.insert(0, str(SOURCE))
spec = importlib.util.spec_from_file_location("phase_bypass", SOURCE / "phase_gate_bypass.py")
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)
D = B.D
M = D.m25
THRESHOLDS = {"region": .5, "start": .5, "stop": .5, "donor": .5, "acceptor": .5}


def example(phase_fail=False, internal_stop=False):
    sequence = list("A"*40)
    sequence[5:12] = list("ATGAAAA")
    sequence[12:14] = list("GT")
    sequence[18:20] = list("AG")
    sequence[20:28] = list("AAAAATAA")
    if internal_stop:
        sequence[8:11] = list("TAA")
    states = np.array([M.I]*5+[M.C]*7+[M.G]*8+[M.C]*8+[M.I]*12)
    region = np.full((40,3), -8., dtype=np.float16)
    region[np.arange(40), states] = 8
    boundary = np.full((40,4), -20., dtype=np.float16)
    for position, channel in ((5,0),(25,1),(12,2),(18,3)):
        boundary[position,channel] = 10
    phase = np.zeros((40,4), dtype=np.float16)
    phase[5,1] = 10
    phase[20,1 if phase_fail else 3] = 10
    return "".join(sequence), region, boundary, phase


def replay(sequence, region, boundary, phase, thresholds=THRESHOLDS, reverse=False):
    _, traces, _, _ = D.trace_decode_orientation(sequence, region, boundary, phase, thresholds, reverse_mapped=reverse)
    assert len(traces) == 1
    t = traces[0]
    probabilities = M._sigmoid(boundary)
    positions = B.P.required_positions(t)
    scores = {p: probabilities[p] for p in positions}
    classes = {p: int(phase[p].argmax()) for p in positions}
    return t, scores, classes


@pytest.mark.parametrize("phase_fail,reverse", [(False,False),(False,True),(True,False),(True,True)])
def test_original_matches_real_frozen_trace_and_production(phase_fail, reverse):
    seq, region, boundary, phase = example(phase_fail)
    trace, scores, classes = replay(seq, region, boundary, phase, reverse=reverse)
    result = B.evaluate_fixed(trace, seq, scores, classes, THRESHOLDS, reverse=reverse)
    B.assert_original(trace, result)


def test_bypass_only_phase_restores_geometry_derived_phase():
    seq, region, boundary, phase = example(True)
    trace, scores, classes = replay(seq, region, boundary, phase)
    result = B.evaluate_fixed(trace, seq, scores, classes, THRESHOLDS, bypass=True)
    assert result["terminal"] == "emitted"
    assert result["model"]["cds"] == [(5,12),(20,28)]
    assert result["model"]["phase"] == [0,2]
    assert result["phase_checks"][-1]["predicted"] == 1
    assert result["phase_checks"][-1]["expected"] == 3


def test_bypass_retains_ORF_stop_and_boundary_thresholds():
    for internal_stop, thresholds, expected in (
        (True, THRESHOLDS, "complete_ORF_internal_stop_check"),
        (False, dict(THRESHOLDS, donor=.99999), "boundary_threshold_filter")):
        seq, region, boundary, phase = example(True, internal_stop)
        trace, scores, classes = replay(seq, region, boundary, phase, thresholds)
        a = B.evaluate_fixed(trace, seq, scores, classes, thresholds)
        B.assert_original(trace, a)
        b = B.evaluate_fixed(trace, seq, scores, classes, thresholds, bypass=True)
        assert b["terminal"] == expected
        assert b["model"] is None


def test_missing_score_and_replay_drift_are_fatal():
    seq, region, boundary, phase = example()
    trace, scores, classes = replay(seq, region, boundary, phase)
    del scores[5]
    with pytest.raises(KeyError):
        B.evaluate_fixed(trace, seq, scores, classes, THRESHOLDS)
    trace, scores, classes = replay(seq, region, boundary, phase)
    result = B.evaluate_fixed(trace, seq, scores, classes, THRESHOLDS)
    changed = deepcopy(trace)
    changed["phase_checks"][0]["predicted"] = 2
    with pytest.raises(AssertionError, match="prefix"):
        B.assert_original(changed, result)


@pytest.mark.parametrize("reverse,chosen_start", [(False,5),(True,11)])
def test_original_orientation_dependent_motif_tie_break(reverse, chosen_start):
    sequence = list("A"*60)
    sequence[5:8] = "ATG"
    sequence[11:14] = "ATG"
    sequence[38:41] = "TAA"
    sequence = "".join(sequence)
    states = np.array([M.I]*8 + [M.C]*33 + [M.I]*19)
    region = np.full((60,3), -8., dtype=np.float16)
    region[np.arange(60),states] = 8
    boundary = np.zeros((60,4), dtype=np.float16)
    phase = np.zeros((60,4), dtype=np.float16)
    phase[:,1] = 10
    trace, scores, classes = replay(sequence, region, boundary, phase, reverse=reverse)
    result = B.evaluate_fixed(trace, sequence, scores, classes, THRESHOLDS, reverse=reverse)
    B.assert_original(trace, result)
    assert result["chosen_positions"]["start"] == [chosen_start]
