from copy import deepcopy
from pathlib import Path
import sys
import unittest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import readout_context_chain_train as D
import verify_context_chain_train as V
from src.m28.chain_set import chain_set_loss
from src.m28.context_chain import ContextChain


class ContextChainReadoutTests(unittest.TestCase):
    def test_saved_components_use_exact_FP32_rounding_and_zero_gain(self):
        value = dict(chains=[[[0,6]]], scores=[0.], additive=[16777216.],
                     exon_additive=[16777216.], junction_additive=[1.], nonadditive=[-16777216.], null=0.)
        V.verify_components(value)
        with self.assertRaises(AssertionError):
            V.verify_components({**value, 'scores':[1.]})
        with self.assertRaises(ValueError):
            V.verify_components({**value, 'nonadditive':[]})
        V.verify_components(dict(chains=[], scores=[], additive=[], exon_additive=[], junction_additive=[],
                                 nonadditive=[], null=.25))

    def test_original_labels_unknown_competition_and_geometry_identity(self):
        species = 'arabidopsis_thaliana'; window = dict(id='fixture', species=species)
        chains = [[[0,12]], [[0,6]], [[18,24]]]; labels = [-1,1,0]
        old = dict(chains=chains, scores=[1.,2.,-1.], selected=[1])
        previous = [dict(window=window, variants={'A31_A3':old})]
        labelled = [dict(window=window, labels={'A31_A3':labels},
                        counts={'A31_A3':D.R.selected_counts(old, labels)},
                        eligible_references=2, eligible_five_or_more=0)]
        new = dict(window=window, chains=chains, scores=[5.,4.,-.1], selected=[0], reference_used_in_scoring=False)
        rows, counts, references = D.paired_readout([new], previous, labelled)
        self.assertEqual(counts[species]['A37']['selected_unknown'], 1)
        self.assertEqual(counts[species]['A37']['selected_known_nonmatch'], 0)
        self.assertEqual(counts[species]['A37']['selected_correct'], 0)
        self.assertEqual(rows[0]['labels'], labels)
        self.assertEqual(references[species]['eligible'], 2)
        for modified in ({**new, 'chains':chains[::-1]}, {**new, 'selected':[1]}, {**new, 'reference_used_in_scoring':True}):
            with self.assertRaises(ValueError):
                D.paired_readout([modified], previous, labelled)

    def test_five_conditions_are_strict_per_species_against_stronger_baseline(self):
        counts = {}
        for species in ('arabidopsis_thaliana', 'oryza_sativa'):
            old = dict(selected_correct=10, selected_known_nonmatch=3, selected_five_or_more_correct=2)
            counts[species] = dict(A31_A3=old, A37={**old, 'selected_correct':11, 'selected_five_or_more_correct':3})
        self.assertEqual(sum(D.gate(counts).values()), 5)
        changed = deepcopy(counts); changed['oryza_sativa']['A37']['selected_correct'] = 10
        self.assertFalse(D.gate(changed)['oryza_sativa_natural_selected_correct_strictly_increased'])
        changed = deepcopy(counts); changed['arabidopsis_thaliana']['A37']['selected_known_nonmatch'] = 4
        self.assertFalse(D.gate(changed)['arabidopsis_thaliana_natural_selected_known_nonmatch_not_increased'])
        changed = deepcopy(counts)
        for values in changed.values():
            values['A37']['selected_five_or_more_correct'] = 2
        self.assertFalse(D.gate(changed)['pooled_natural_selected_five_or_more_correct_strictly_increased'])

    def test_branch_exposure_matches_actual_toy_optimizer_including_unknown_calls(self):
        torch.manual_seed(0); model = ContextChain(8,16)
        features, dna = torch.randn(1,5,8), torch.randn(1,30,5)
        cases = (([],[]), ([[[0,6]]],[-1]), ([[[0,6],[12,18]]],[-1]),
                 ([[[0,6]], [[0,6],[12,18]]],[1,-1]))
        records = {str(i):dict(candidates=[dict(chain=c, label=y) for c, y in zip(chains, labels)])
                   for i, (chains, labels) in enumerate(cases)}
        draws = [dict(window_id=str(i)) for i in range(4)]
        optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01)
        for _ in range(2):
            for chains, labels in cases:
                optimizer.zero_grad(set_to_none=True)
                score = model(features, dna, [30], chains)
                chain_set_loss(score, chains, labels)['loss'].backward(); optimizer.step()
        actual = {name:int(optimizer.state[p]['step']) if p in optimizer.state else 0
                  for name, p in model.named_parameters()}
        expected = V.expected_optimizer_steps(list(actual), draws, records, passes=2)
        self.assertEqual(actual, expected)
        self.assertEqual(actual['head.null.weight'], 8)
        self.assertEqual(actual['head.exon.0.weight'], 6)
        self.assertEqual(actual['head.junction.0.weight'], 4)


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
