"""compact offset == parity-gate build_logits (A_full), bit for bit, on a real 400-kb TRAIN region, both strands."""
import sys
import numpy as np
import pytest
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_gate_parity as GP                 # noqa: E402
from constants import MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR   # noqa: E402
import a1_offset as O                        # noqa: E402

SEQ = 'NC_003070.9'


@pytest.mark.parametrize('strand', '+-')
def test_compact_equals_gate(strand):
    rid, sequence, pos, neg = GP.load_tracks(SEQ)
    a, b = 6_000_000, 6_400_000
    if strand == '+':
        s, p = sequence[a:b], pos[a:b]
        tt = np.load(K.CACHE / f'tt_{SEQ}_+.npy')[a:b]
    else:
        L = len(sequence)
        s = O.rc(sequence)[L - b:L - a]
        p = np.ascontiguousarray(neg[::-1])[L - b:L - a]
        tt = np.load(K.CACHE / f'tt_{SEQ}_-.npy')[L - b:L - a]
    ref = GP.build_logits(s, p, tt, 'A_full')
    got = O.expand(O.compact(s, p, tt, MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR))
    assert got.dtype == ref.dtype == np.float32
    assert np.array_equal(got, ref)
    assert np.count_nonzero(ref[:, 13]) > 0 and np.count_nonzero(ref[:, 7]) > 0


def test_context_slice_minus_orientation():
    rng = np.random.default_rng(0)
    L = 5000
    c = {'+': rng.standard_normal((5, L)).astype(np.float32), '-': rng.standard_normal((5, L)).astype(np.float32)}
    a, b = 1200, 2100
    assert np.array_equal(O.context_offset(c, '+', a, b, L), O.expand(c['+'], a, b))
    assert np.array_equal(O.context_offset(c, '-', a, b, L), O.expand(c['-'], L - b, L - a))
