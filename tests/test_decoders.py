"""Correctness checks for TA-DECODER-M3 structured decoders (run where torch is available:
coding-rna env on baobab). Verifies CRF nll/viterbi + constrained_decode fragmentation fix."""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from src.screen_anchor.decoders import LinearChainCRF, LinearChainCRFVec, constrained_decode
from src.screen_anchor.data import CLASS_CDS, CLASS_GENEBODY_NC


def _copy_params(ref, vec):
    with torch.no_grad():
        vec.trans.copy_(ref.trans); vec.start.copy_(ref.start); vec.end.copy_(ref.end)


def test_crfvec_matches_reference_partition_and_gold():
    """Vectorized (parallel-scan) CRF must equal the reference O(W) CRF on partition + gold +
    nll, for full AND padded batches (this is the goal's CK1 'vectorized==reference' gate)."""
    torch.manual_seed(0)
    for W, pad in [(40, 0), (50, 7), (1, 0)]:
        B, C = 4, 3
        em = torch.randn(B, W, C)
        Y = torch.randint(0, 3, (B, W))
        mask = torch.ones(B, W, dtype=torch.bool)
        if pad:
            Y[:, -pad:] = -100
            mask[:, -pad:] = False
        ref, vec = LinearChainCRF(C), LinearChainCRFVec(C)
        _copy_params(ref, vec)
        assert torch.allclose(ref._forward_alg(em, mask), vec._forward_alg(em, mask), atol=1e-4), f"partition W={W}"
        assert torch.allclose(ref._gold_score(em, Y, mask), vec._gold_score(em, Y, mask), atol=1e-4), f"gold W={W}"
        assert torch.allclose(ref.nll(em, Y, mask), vec.nll(em, Y, mask), atol=1e-4), f"nll W={W}"


def test_crfvec_partition_ge_gold_and_backprop():
    torch.manual_seed(3)
    B, W, C = 4, 64, 3
    em = torch.randn(B, W, C, requires_grad=True)
    Y = torch.randint(0, 3, (B, W)); mask = torch.ones(B, W, dtype=torch.bool)
    crf = LinearChainCRFVec(C)
    nll = crf.nll(em, Y, mask)
    assert torch.isfinite(nll); assert (crf._forward_alg(em.detach(), mask) + 1e-4 >= crf._gold_score(em.detach(), Y, mask)).all()
    nll.backward(); assert em.grad is not None and torch.isfinite(em.grad).all()


def test_crf_nll_finite_and_backprop():
    torch.manual_seed(0)
    B, W, C = 4, 50, 3
    em = torch.randn(B, W, C, requires_grad=True)
    Y = torch.randint(0, 3, (B, W))
    Y[:, -5:] = -100                      # padding
    mask = Y != -100
    crf = LinearChainCRF(C)
    nll = crf.nll(em, Y, mask)
    assert torch.isfinite(nll)
    nll.backward()
    assert em.grad is not None and torch.isfinite(em.grad).all()


def test_crf_partition_ge_gold():
    # logZ (marginal over all paths) must be >= any single gold path score
    torch.manual_seed(1)
    B, W, C = 3, 40, 3
    em = torch.randn(B, W, C)
    Y = torch.randint(0, 3, (B, W))
    mask = torch.ones(B, W, dtype=torch.bool)
    crf = LinearChainCRF(C)
    logZ = crf._forward_alg(em, mask)
    gold = crf._gold_score(em, Y, mask)
    assert (logZ + 1e-4 >= gold).all()


def test_crf_viterbi_shape_and_range():
    torch.manual_seed(2)
    B, W, C = 2, 30, 3
    em = torch.randn(B, W, C)
    mask = torch.ones(B, W, dtype=torch.bool)
    path = LinearChainCRF(C).viterbi(em, mask)
    assert path.shape == (B, W)
    assert int(path.min()) >= 0 and int(path.max()) <= C - 1


def test_constrained_decode_reduces_fragmentation():
    # a CDS run broken by one tiny intergenic gap -> should merge into one gene-body
    arr = np.ones(60, dtype=np.int8) * CLASS_CDS
    arr[30] = 0                               # 1-base intergenic gap splits the CDS
    before_genes = int((np.diff((arr > 0).astype(int)) > 0).sum()) + int(arr[0] > 0)
    out = constrained_decode({"s": arr}, min_cds_len=5, max_fill_gap=20)["s"]
    after_genes = int((np.diff((out > 0).astype(int)) > 0).sum()) + int(out[0] > 0)
    assert after_genes <= before_genes        # gap filled -> fewer/equal genes
    assert (out > 0).all()                     # the whole region is now one gene-body


def test_constrained_decode_drops_tiny_cds():
    arr = np.zeros(100, dtype=np.int8)
    arr[10:13] = CLASS_CDS                     # 3-base CDS fragment (< min_cds_len)
    out = constrained_decode({"s": arr}, min_cds_len=30, max_fill_gap=0)["s"]
    assert (out == CLASS_CDS).sum() == 0       # tiny CDS dropped
    assert (out[10:13] == CLASS_GENEBODY_NC).all()
