"""Artificial paired TRAIN-score diagnostics, no actual weights or inputs."""
from copy import deepcopy
import unittest
import torch

import diagnose_paired_emission_mechanism as M
from test_m28_local_structured_raw import receipts


class PairedMechanismTests(unittest.TestCase):
    def test_numeric_sources_allow_scientific_FAIL_but_not_stale_or_unverified_results(self):
        _, prepared, checked = receipts()
        M.validate_sources(prepared, checked)
        for report in (prepared, checked):
            for values in report['selected_counts_per_species'].values():
                values['A44_JOINT'] = deepcopy(values['A44_LOCAL'])
            report['predeclared_TRAIN_gates'], report['TRAIN_preparation_supported'] = M.D.gates(report['selected_counts_per_species'])
        self.assertFalse(prepared['TRAIN_preparation_supported'])
        M.validate_sources(prepared, checked)
        for change in ('old_source', 'scope', 'proof', 'fake_PASS'):
            a, b = deepcopy(prepared), deepcopy(checked)
            if change == 'old_source': b['source_jobs']['readout'] = '13424410'
            elif change == 'scope': a['DEV_used'] = True
            elif change == 'proof': b[M.I.CPU_PROOFS[0]] = False
            else:
                a['TRAIN_preparation_supported'] = b['TRAIN_preparation_supported'] = True
            with self.subTest(change=change), self.assertRaises(ValueError): M.validate_sources(a, b)

    def test_literal_six_role_terms_conserve_FP64_gain_and_shared_logit_shift_invariance(self):
        chain = [[1, 7], [10, 16]]; local = torch.zeros(18, 15)
        # Literal original start/CDS/donor/intron/acceptor/stop precedence.
        roles = {1: 13, 2: 13, 3: 13, 4: 1, 5: 2, 6: 9,
                 7: 6, 8: 6, 9: 6, 10: 10, 11: 2, 12: 3, 13: 14, 14: 14, 15: 14}
        for position, role in roles.items(): local[position, role] = role
        terms = M.components(local, chain)
        self.assertEqual(terms, dict(coding_phase=8., intron_phase=18., donor_phase=9.,
                                    acceptor_phase=10., start=39., stop=42.))
        self.assertEqual(sum(terms.values()), 126.)
        offsets = torch.arange(18, dtype=torch.float32)[:, None]*3
        self.assertEqual(M.components(local+offsets, chain), terms)

    def test_unknown_and_nonoverlapping_nonmatches_do_not_manufacture_conflict_margin(self):
        gold, bad, unknown, distant = [[1, 7]], [[1, 10]], [[1, 13]], [[16, 22]]
        candidates = [dict(chain=c, label=y) for c, y in ((gold, 1), (bad, 0), (unknown, -1), (distant, 0))]
        local = torch.zeros(24, 15); joint = torch.ones(24, 15)
        local[:, 13] = 2.; local[:, 14] = 3.; joint[:, 13] = 4.; joint[:, 14] = 5.
        local[7:13, 1:4] = 1e5; joint[7:13, 1:4] = 1e5
        predictions = {a: dict(chains=[]) for a in M.ARMS}
        row = M.diagnose_window(candidates, dict(zip(M.ARMS, (local, joint))), predictions)[0]
        self.assertEqual(row['finite_known_conflict_count'], 1)
        self.assertEqual({row['arms'][a]['finite_best_known_conflict_index'] for a in M.ARMS}, {1})
        self.assertEqual(row['selection_transition'], 'neither_selected')
        self.assertTrue(all(row['arms'][a]['gain_positive'] for a in M.ARMS))
        no_conflict = M.diagnose_window([candidates[0], candidates[2], candidates[3]],
            dict(zip(M.ARMS, (local, joint))), predictions)[0]
        self.assertIsNone(no_conflict['joint_minus_local_finite_known_conflict_margin'])
        self.assertIsNone(no_conflict['arms'][M.ARMS[0]]['finite_best_known_conflict_index'])

    def test_exact_retention_transitions_do_not_depend_on_gain_or_alternative_chain(self):
        chains = [[[i, i+6]] for i in (1, 10, 19, 28)]
        candidates = [dict(chain=c, label=1) for c in chains]
        potentials = {a: torch.zeros(40, 15) for a in M.ARMS}
        for local in potentials.values(): local[:, 13:15] = 1.
        predictions = {M.ARMS[0]: dict(chains=[chains[0], chains[2]]),
                       M.ARMS[1]: dict(chains=[chains[0], chains[1], [[28, 37]]])}
        rows = M.diagnose_window(candidates, potentials, predictions)
        self.assertEqual([r['selection_transition'] for r in rows],
                         ['both_selected', 'newly_selected', 'lost_selected', 'neither_selected'])
        summary = M.summarize(rows)
        self.assertEqual(summary['reference_pairs'], 4)
        self.assertEqual(summary['selection_transitions'], dict(both_selected=1, newly_selected=1, lost_selected=1, neither_selected=1))
        self.assertEqual(summary['references_without_finite_known_conflict'], 4)
        self.assertEqual(summary['joint_minus_local_finite_known_conflict_margin']['n'], 0)
        self.assertIsNone(summary['joint_minus_local_finite_known_conflict_margin']['mean'])
        for arm in M.ARMS:
            self.assertEqual(summary['arms'][arm]['selected'], 2)
            self.assertEqual(summary['arms'][arm]['positive_not_selected'], 2)
            self.assertEqual(summary['arms'][arm]['nonpositive_not_selected'], 0)


if __name__ == '__main__': unittest.main()
