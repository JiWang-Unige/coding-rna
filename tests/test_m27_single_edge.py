import sys
from pathlib import Path
import numpy as np
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts/experiments/M27-ALLELIC-INTEGRITY"))
sys.path.insert(0, str(ROOT/"refs/repos/annevo-2026"))
import edge_decode as E
from src import HMM as H


def toy():
    m = np.full((5, 3, 3), -np.inf, dtype=np.float32)
    m[:, 0, :] = 0
    m[:, 1, 0] = 0
    m[:, 2, 0] = 0
    return m, np.zeros(4, dtype=np.uint8), {"CDS1": 0, "CDS1_TA": 1, "CDS2": 2}


def test_exact_one_entry_one_position_and_no_mutation():
    m, c, s = toy()
    old_m, old_c = m.copy(), c.copy()
    new, codes = E.release_edge(m, c, 2, s)
    assert np.array_equal(m, old_m) and np.array_equal(c, old_c)
    assert np.array_equal(new[:5], m)
    assert np.argwhere(new[5] != m[0]).tolist() == [[1, 2]]
    assert np.flatnonzero(codes != c).tolist() == [2]
    assert new[5, 1, 2] == m[0, 0, 2] == 0


def test_actual_native_numba_path_uses_only_released_edge():
    m, c, s = toy()
    emissions = np.full((4, 3), -10, dtype=np.float32)
    emissions[1, 1] = emissions[2, 2] = emissions[3, 0] = 0
    baseline = H._viterbi_core_numba_float32(emissions, m, c)
    new, codes = E.release_edge(m, c, 2, s)
    relaxed = H._viterbi_core_numba_float32(emissions, new, codes)
    assert relaxed.tolist() == [0, 1, 2, 0]
    assert baseline.tolist() != relaxed.tolist()


def test_native_annevo_a_matrix_only_cds_continuation_changes():
    s, n = H.define_state(1)
    init = np.full((n, n), -np.inf, dtype=np.float32)
    init = H.set_transition_matrix_common_state(init, s, 1, 0, 0, 0, 0)
    d = H.set_transition_matrix_conditional_state(init, s, 1, 0, 0, 0, 0)
    m = np.stack([d[b] for b in "ATCGN"])
    new, _ = E.release_edge(m, np.zeros(4, dtype=np.uint8), 2, s)
    assert np.argwhere(new[-1] != m[0]).tolist() == [[s["CDS1_TA"], s["CDS2"]]]
    for source, dest in (("ASS1_TA", "CDS2"), ("CDS1_TA", "DSS2"), ("ASS1_TA", "DSS2")):
        assert np.isneginf(new[-1, s[source], s[dest]])


def test_wrong_base_or_already_allowed_edge_stops():
    m, c, s = toy()
    c[2] = 1
    with unittest.TestCase().assertRaises(RuntimeError):
        E.release_edge(m, c, 2, s)
    c[2] = 0
    m[0, 1, 2] = 0
    with unittest.TestCase().assertRaises(RuntimeError):
        E.release_edge(m, c, 2, s)


def test_target_coordinate_and_strand_scope():
    region = (E.SPAN[0]-100, E.SPAN[1]+100, "chr22", 1)
    assert E.target_index(region) == E.POS-region[0]
    assert E.target_index((*region[:3], -1)) is None
    assert E.target_index((*region[:2], "chr21", 1)) is None
    with unittest.TestCase().assertRaisesRegex(RuntimeError, "anchors"):
        E.target_index((E.POS-5, E.POS+5, "chr22", 1))


def test_component_sum_includes_noncoding_and_transitions_but_not_initial_emission():
    e = np.array([[999, 999], [2, 20], [3, 30], [4, 40]], dtype=np.float32)
    cols, path, codes = np.array([0, 1]), np.array([0, 1, 0, 0]), np.zeros(4, dtype=np.uint8)
    m = np.ones((5, 2, 2), dtype=np.float32)
    score = E.components(e, cols, path, m, codes)
    assert score == {"positions": 3, "emission": 27.0, "transition": 3.0, "total": 30.0,
                     "impossible_transition_indices": []}
    m[0, 1, 0] = -np.inf
    score = E.components(e, cols, path, m, codes)
    assert score["total"] is None and score["impossible_transition_indices"] == [2]


if __name__ == "__main__":
    suite = unittest.TestSuite(unittest.FunctionTestCase(v) for k, v in sorted(globals().copy().items()) if k.startswith("test_"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
