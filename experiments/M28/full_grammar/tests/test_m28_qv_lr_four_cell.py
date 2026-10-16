"""Artificial accounting checks, no biological input or model dependencies."""
import copy
import unittest
import compare_full_grammar_qv_lr_four_cell as F


def prf(tp, predicted, reference):
    return dict(tp=tp, predicted=predicted, reference=reference, precision=tp/predicted,
                recall=tp/reference, f1=2*tp/(predicted+reference))


def fixture():
    refs = (6000, 6000, 5739, 1600, 1600, 1749)  # Artificial partition, original fixed total only.
    rates = {'lora_low':(5, 8), 'lora_high':(7, 6), 'dense_low':(4, 7), 'dense_high':(7, 5)}
    reports = {}
    for index, (name, factors) in enumerate(rates.items()):
        rows = {}
        for (species, seqid, length), reference in zip(F.SCOPE, refs):
            tp = reference*factors[F.SPECIES.index(species)]//10
            bg = dict(bp=100000, predicted_span_bp=(index+1)*1000, predicted_CDS_bp=(index+1)*500,
                      wholly_background_unique_chains=index+1, span_FPR=(index+1)/100.,
                      wholly_background_chains_per_Mb=(index+1)*10.)
            rows[seqid] = dict(species=species, seqid=seqid, length=length,
                reference_primary_records=reference, reference_assessable_isoform_records=reference+2,
                retained_unsupported_reference_strata={'longer_than_native_window':1},
                exact_chain=prf(tp, reference, reference), CDS_base_unstranded_union=prf(800, 1000, 1000),
                background=bg, prediction_records=reference, duplicate_chain_records=0,
                exact_other_assessable_isoform_not_primary_TP=1, phase_checked_prediction_records=reference,
                prediction_nonexact_primary=reference-tp)
        def aggregate(parts):
            bg = {k:sum(r['background'][k] for r in parts) for k in F.BG_COUNTS}
            bg.update(span_FPR=bg['predicted_span_bp']/bg['bp'],
                      wholly_background_chains_per_Mb=bg['wholly_background_unique_chains']*1e6/bg['bp'])
            return dict(**{e:prf(*(sum(r[e][k] for r in parts) for k in ('tp', 'predicted', 'reference'))) for e in F.ENDPOINTS},
                background=bg, **{k:sum(r[k] for r in parts) for k in F.RECORD_COUNTS})
        per_species = {s:aggregate([r for r in rows.values() if r['species'] == s]) for s in F.SPECIES}
        f1s = [per_species[s]['exact_chain']['f1'] for s in F.SPECIES]
        reports[name] = dict(reference_policy='artificial_common_reference', chromosomes=rows,
            per_species=per_species, pooled=aggregate(list(rows.values())),
            macro_species_exact_chain_F1=sum(f1s)/2, worst_species_exact_chain_F1=min(f1s))
    return reports


class FourCellAccountingTests(unittest.TestCase):
    def test_signed_interaction_and_per_cell_worst_before_difference(self):
        reports = fixture()
        scopes = F.summarize(reports)
        key = 'exact_chain.f1_percent'
        worst = scopes['worst_species']['contrasts']['lora_high_minus_low'][key]
        per_species = [scopes[s]['contrasts']['lora_high_minus_low'][key] for s in F.SPECIES]
        self.assertGreater(worst, 9.)
        self.assertLess(min(per_species), -19.)
        pooled = scopes['pooled']
        effects = pooled['contrasts']
        self.assertAlmostEqual(effects['difference_in_differences'][key],
            effects['dense_high_minus_low'][key]-effects['lora_high_minus_low'][key])
        tp = {name:report['pooled']['exact_chain']['tp'] for name, report in reports.items()}
        self.assertAlmostEqual(effects['difference_in_differences'][key],
            100*(tp['dense_high']-tp['dense_low']-tp['lora_high']+tp['lora_low'])/22688)
        self.assertGreater(effects['difference_in_differences'][key], 7.)
        self.assertAlmostEqual(effects['lora_high_minus_low']['background.span_FPR_percent'], 1.)
        self.assertEqual(pooled['arms']['lora_low']['phase_checked_prediction_records'], 22688)
        self.assertEqual(pooled['arms']['lora_low']['nonprimary_unique_chains'],
            22688-reports['lora_low']['pooled']['exact_chain']['tp'])

    def test_reject_missing_cell_changed_reference_geometry_and_bad_prf(self):
        reports = fixture()
        missing = copy.deepcopy(reports); del missing['dense_high']
        with self.assertRaises(ValueError):
            F.summarize(missing)
        for key, value in (('length', 12), ('retained_unsupported_reference_strata', {})):
            wrong = copy.deepcopy(reports)
            wrong['dense_high']['chromosomes']['NC_003070.9'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                F.summarize(wrong)
        wrong = copy.deepcopy(reports)
        wrong['dense_high']['per_species'][F.SPECIES[0]]['exact_chain']['f1'] += .01
        with self.assertRaises(ValueError):
            F.summarize(wrong)

    def test_reject_background_phase_nonprimary_and_pooled_count_mismatch(self):
        reports = fixture()
        for kind in ('background', 'phase', 'nonprimary', 'pooled'):
            wrong = copy.deepcopy(reports)
            row = wrong['dense_high']['chromosomes']['NC_003070.9']
            if kind == 'background':
                row['background']['span_FPR'] += .01
            elif kind == 'phase':
                row['phase_checked_prediction_records'] -= 1
            elif kind == 'nonprimary':
                row['prediction_nonexact_primary'] -= 1
            else:
                wrong['dense_high']['pooled']['background']['wholly_background_unique_chains'] += 1
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                F.summarize(wrong)
