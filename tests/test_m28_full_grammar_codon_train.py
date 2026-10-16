"""Artificial identity failures and paired count arithmetic; no research data."""
import copy
import unittest
import torch
import infer_full_grammar_codon_train as D
from compare_full_grammar_codon_train import paired_scopes, SPECIES, ENDPOINTS
from compare_full_grammar_uniform_coarse_train import checked_prf


def prf(tp, predicted, reference):
    return dict(tp=tp, predicted=predicted, reference=reference,
        precision=tp/predicted if predicted else 0., recall=tp/reference if reference else 0.,
        f1=2*tp/(predicted+reference) if predicted+reference else 0.)


def report_fixture(scores):
    return dict(inference_complete=True, checkpoint=str(D.CHECKPOINT), lengths={'toy14':14},
        reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
        optimizer_steps_this_run=0, backend=D.BACKEND, score_kind=D.SCORE_KIND,
        raw_codon_scores=scores.tolist(), centered_codon_scores=D.centered_codon_scores(scores).tolist(),
        configuration=dict(D.A.IDENTITY), scopes=[dict(seqid='toy14', strand=strand, windows=[[0, 14]],
            reference_used=False, joint_whole_oriented_contig_decode=True, context_policy=D.I.R.CONTEXT,
            pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
            score_accumulation=D.ACCUMULATION, backend=D.BACKEND) for strand in ('+', '-')])


def paired_fixture():
    reports = {}
    references = (3000, 3000, 3000, 4504, 4504, 4680)
    for name in ('baseline', 'codon'):
        chromosomes = {}
        for (species, seqid, length), n in zip(D.I.SCOPE, references):
            fraction = (7 if name == 'baseline' else 6) if species == SPECIES[0] else (4 if name == 'baseline' else 5)
            tp = n*fraction//8
            chromosomes[seqid] = dict(species=species, seqid=seqid, length=length,
                reference_primary_records=n, reference_assessable_isoform_records=n+1,
                retained_unsupported_reference_strata={'unsupported_fixture':1}, background={},
                exact_chain=prf(tp, n, n), CDS_base_unstranded_union=prf(tp*3, n*3, n*3))
        def aggregate(rows):
            return dict(background={}, **{e:prf(*(sum(r[e][k] for r in rows) for k in ('tp', 'predicted', 'reference'))) for e in ENDPOINTS})
        reports[name] = dict(reference_policy='fixed_fixture_reference_policy', chromosomes=chromosomes,
            per_species={s:aggregate([r for r in chromosomes.values() if r['species'] == s]) for s in SPECIES},
            pooled=aggregate(list(chromosomes.values())))
    return reports


class CodonTrainEndpointTests(unittest.TestCase):
    def test_raw_requires_actual_codon_backend_and_both_original_strands(self):
        scores = torch.arange(64, dtype=torch.float32)/100
        good = report_fixture(scores)
        D.validate_raw_report(good, 'toy14', 14, scores)
        for key, value in (('backend', 'exact_exon_budget_skyline'), ('score_kind', 'role_only'),
                           ('checkpoint', 'wrong_checkpoint'), ('reference_used', True)):
            changed = copy.deepcopy(good)
            changed[key] = value
            with self.assertRaises(ValueError):
                D.validate_raw_report(changed, 'toy14', 14, scores)
        changed = copy.deepcopy(good)
        changed['scopes'][1]['strand'] = '+'
        with self.assertRaises(ValueError):
            D.validate_raw_report(changed, 'toy14', 14, scores)

    def test_coefficients_centering_order_and_coverage_mismatch_are_rejected(self):
        scores = torch.arange(64, dtype=torch.float32)/100
        good = report_fixture(scores)
        for key in ('raw_codon_scores', 'centered_codon_scores'):
            changed = copy.deepcopy(good)
            changed[key][14] += .1
            with self.assertRaises(ValueError):
                D.validate_raw_report(changed, 'toy14', 14, scores)
        for key, value in (('windows', [[0, 13]]), ('score_accumulation', 'FP64_per_base_forward_order')):
            changed = copy.deepcopy(good)
            changed['scopes'][0][key] = value
            with self.assertRaises(ValueError):
                D.validate_raw_report(changed, 'toy14', 14, scores)
        changed = copy.deepcopy(good)
        changed['configuration']['codon_order'] = list(reversed(D.A.CODONS))
        with self.assertRaises(ValueError):
            D.validate_raw_report(changed, 'toy14', 14, scores)

    def test_counts_reconstruct_and_worst_species_is_not_worst_paired_difference(self):
        scopes = paired_scopes(paired_fixture())
        self.assertEqual(scopes['pooled']['arms']['codon']['exact_chain']['reference'], 22688)
        self.assertEqual(scopes['macro_species']['codon_minus_baseline_PRF_percentage_points']['exact_chain']['f1'], 0.)
        self.assertEqual(scopes['worst_species']['codon_minus_baseline_PRF_percentage_points']['exact_chain']['f1'], 12.5)
        worst_difference = min(scopes[s]['codon_minus_baseline']['exact_chain']['PRF_percentage_points']['f1'] for s in SPECIES)
        self.assertEqual(worst_difference, -12.5)
        self.assertNotEqual(worst_difference, scopes['worst_species']['codon_minus_baseline_PRF_percentage_points']['exact_chain']['f1'])

    def test_changed_reference_population_or_aggregate_counts_are_rejected(self):
        reports = paired_fixture()
        seqid = D.I.SCOPE[0][1]
        changed = copy.deepcopy(reports)
        changed['codon']['chromosomes'][seqid]['reference_primary_records'] -= 1
        with self.assertRaises(ValueError):
            paired_scopes(changed)
        changed = copy.deepcopy(reports)
        changed['codon']['pooled']['exact_chain'] = prf(1, 22688, 22688)
        with self.assertRaises(ValueError):
            paired_scopes(changed)
        invalid = prf(5, 10, 10)
        invalid['f1'] = .8
        with self.assertRaises(ValueError):
            checked_prf(invalid)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
