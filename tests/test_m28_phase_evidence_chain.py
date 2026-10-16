"""A26 phase, orientation, evidence flow and unchanged unknown semantics."""
import math
import unittest
import torch
from src.m28.labels import _render
from src.m28.junction_chain import JunctionChainHead
from src.m28.phase_evidence_chain import candidate_roles, phase_evidence, PhaseEvidenceChainHead
from src.m28.chain_set import chain_set_loss


class PhaseEvidenceTests(unittest.TestCase):
    def test_phase_carry_boundary_precedence_and_both_orientations(self):
        chain = ((2, 6), (11, 22))
        actual = candidate_roles(chain, 24)
        expected = [0, 0, 13, 13, 13, 7]+[4]*5+[11, 3, 1, 2, 3, 1, 2, 3]+[14]*3+[0, 0]
        self.assertEqual(actual.tolist(), expected)
        for strand, cds in (('+', [(102, 106, 0), (111, 122, 2)]),
                            ('-', [(102, 113, 2), (118, 122, 0)])):
            rendered = _render(dict(CDS=cds, CDS_complete_supported=True),
                               dict(start=100, end=124, width=24, valid_bases=24, strand=strand))
            self.assertEqual(actual.tolist(), rendered['seg'].tolist())

    def test_exact_state_odds_counts_padding_and_per_base_shift_invariance(self):
        chain = ((2, 6), (11, 22))
        roles = candidate_roles(chain, 24)
        logits = torch.zeros(24, 15, dtype=torch.float64)
        logits[:, 0] = 2
        for state in range(1, 15):
            logits[roles == state, state] = state+2
        evidence = phase_evidence(logits, [chain])
        for state in range(1, 15):
            count = int((roles == state).sum())
            self.assertEqual(float(evidence[0, state-1]), float(state if count else 0))
            self.assertEqual(float(evidence[0, 14+state-1]), math.log1p(count))
        shifted = logits+torch.arange(24, dtype=torch.float64)[:, None]*10
        torch.testing.assert_close(evidence, phase_evidence(shifted, [chain]), rtol=0, atol=0)
        padded = torch.cat([logits, torch.full((5, 15), 1e9, dtype=torch.float64)])
        torch.testing.assert_close(evidence, phase_evidence(padded, [chain], 24), rtol=0, atol=0)
        logits[[0, 1, 22, 23]] = -1e9
        torch.testing.assert_close(evidence, phase_evidence(logits, [chain]), rtol=0, atol=0)

    def test_zero_prior_preserves_all_fresh_A21_common_parameters_and_scores(self):
        torch.manual_seed(0); old = JunctionChainHead(8).double()
        torch.manual_seed(0); head = PhaseEvidenceChainHead(8).double()
        for key, value in old.state_dict().items():
            self.assertTrue(torch.equal(value, head.state_dict()[key]), key)
        self.assertEqual(sum(p.numel() for p in head.parameters())-sum(p.numel() for p in old.parameters()), 28)
        hidden = torch.randn(24, 8, dtype=torch.float64)
        chains = [((2, 6), (11, 22)), ((2, 20),)]
        prior = old(hidden, chains)
        actual = head(hidden, chains, phase_evidence(torch.randn(24, 15, dtype=torch.float64), chains))
        for key, value in prior.items():
            torch.testing.assert_close(value, actual[key], rtol=0, atol=0)

    def test_candidate_role_perturbation_enters_score_with_frozen_local_logits(self):
        head = PhaseEvidenceChainHead(8).double()
        with torch.no_grad():
            head.path_evidence.weight[0, 3] = 1  # mean intron0 vs intergenic odds
        hidden = torch.randn(24, 8, dtype=torch.float64)
        chains = [((2, 6), (11, 22))]
        logits = torch.randn(24, 15, dtype=torch.float64, requires_grad=True)
        before = head(hidden, chains, phase_evidence(logits, chains))
        changed = logits.detach().clone(); changed[8, 4] += 5
        after = head(hidden, chains, phase_evidence(changed, chains))
        self.assertAlmostEqual(float(after['logits'][0]-before['logits'][0]), 1., places=12)
        before['logits'].sum().backward()
        self.assertIsNone(logits.grad)
        self.assertGreater(float(head.path_evidence.weight.grad.abs().sum()), 0.)

    def test_empty_touching_chains_and_unknown_direct_gradient(self):
        head = PhaseEvidenceChainHead(8).double()
        hidden = torch.randn(24, 8, dtype=torch.float64)
        logits = torch.randn(24, 15, dtype=torch.float64)
        empty = head(hidden, [], phase_evidence(logits, []))
        self.assertEqual(tuple(empty['phase_evidence_term'].shape), (0,))
        loss = chain_set_loss(empty, [], [])['loss']; loss.backward()
        self.assertEqual(float(loss), 0.)
        self.assertTrue((head.path_evidence.weight.grad == 0).all())
        head.zero_grad(set_to_none=True)
        chains = [((2, 6), (11, 22)), ((2, 6), (6, 21))]
        scores = head(hidden, chains, phase_evidence(logits, chains))
        self.assertTrue(torch.isfinite(scores['logits']).all())
        scores['logits'].retain_grad()
        chain_set_loss(scores, chains, [1, -1])['loss'].backward()
        self.assertEqual(float(scores['logits'].grad[1]), 0.)
        self.assertNotEqual(float(scores['logits'].grad[0]), 0.)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
