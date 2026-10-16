from collections import Counter
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_junction_chain as V
from test_m28_junction_chain_raw import receipts
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def decoded_receipt(fitted):
    scopes = []
    for ordinal, (species, lengths) in enumerate(V.I.C.I.DEV_LENGTHS.items()):
        for seqid, length in lengths.items():
            for strand in ('+', '-'):
                scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
                    windows=[[a, min(a+V.I.M.WIDTH, length)] for a in V.I.M.starts(length)],
                    file='contig_'+str(ordinal)+'_'+strand+'.jsonl'))
    return dict(arm='A21', completed=True, inference_complete=True, configuration=V.I.J.CONFIG,
        protocol=V.I.DEV_PROTOCOL, TRAIN_protocol=V.I.J.PROTOCOL, TRAIN_preparation_gate_passed=True,
        checkpoint=fitted['final_checkpoint'], optimizer_steps=4608, generator_checkpoint=V.I.C.GENERATOR,
        feature_revision=V.I.C.I.F.REVISION, feature_spec=V.I.C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA', reference_used=False, labels_or_feature_cache_used=False,
        optimizer_steps_this_run=0, DEV_used_for_selection=False, test_setaria_used=False,
        windows=8690, scopes=scopes,
        lengths={q: n for qs in V.I.C.I.DEV_LENGTHS.values() for q, n in qs.items()})


def raw_row():
    return dict(window_index=0, window_id='synthetic|chr|0|+', chains=[[[0, 9]]], score=[1.],
        score_components=dict(exon_additive=[2.], junction_additive=[0.], nonadditive=[1.]), window_null=2.,
        decoder_sources=[['beam4', 'greedy']], score_kind=V.I.SCORE_KIND,
        quality_logit_threshold=0., reference_used=False,
        counts=dict(free_candidates=1, nonpositive_candidates=0))


class JunctionEvaluationTests(unittest.TestCase):
    def test_full_8690_grid_required_not_explicit_subset_or_different_encoder_context(self):
        fitted, _, _ = receipts(); good = decoded_receipt(fitted); V.validate_complete(fitted, good)
        for key, value in (('scopes', good['scopes'][:-1]), ('input_mode', 'explicit_raw_FASTA'),
                           ('feature_spec', {'encoder_context': 'whole_4096'}), ('reference_used', True)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                V.validate_complete(fitted, {**good, key: value})
        truncated = copy.deepcopy(good); truncated['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError):
            V.validate_complete(fitted, truncated)

    def test_matched_pool_and_actual_components_not_just_plausible_score_receipt(self):
        scope = dict(species='synthetic', seqid='chr', strand='+', windows=[[0, 9]])
        row = raw_row(); old = {**row, 'score': [-1.], 'score_kind': 'A6_complete_chain_logit_minus_window_null'}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); current = root/'new.jsonl'; previous = root/'old.jsonl'
            previous.write_text(json.dumps(old)+'\n'); current.write_text(json.dumps(row)+'\n')
            counts = Counter()
            self.assertEqual(len(list(V.matched_records(current, scope, previous, scope, counts))), 1)
            self.assertEqual(dict(counts), dict(free_candidates=1, nonpositive_candidates=0))
            for key, value in (('chains', [[[0, 6]]]), ('decoder_sources', [['greedy']]),
                               ('window_null', 3.), ('score', [float('nan')]),
                               ('score_components', dict(exon_additive=[2.], junction_additive=[], nonadditive=[1.]))):
                current.write_text(json.dumps({**row, key: value})+'\n')
                with self.subTest(key=key), self.assertRaises(ValueError):
                    list(V.matched_records(current, scope, previous, scope, Counter()))
            current.write_text('')
            with self.assertRaises(ValueError):
                list(V.matched_records(current, scope, previous, scope, Counter()))

    def test_reconstruction_closes_gff_count_and_phase_without_reference_access(self):
        row = raw_row()
        scope = dict(species='synthetic', seqid='chr', strand='+', ordinal=0, windows=[[0, 9]], file='scope.jsonl',
                     quality_counts=row['counts'], assembly_counts=assemble_predictions([row], [[0, 9]], '+')['counts'])
        decoded = dict(scopes=[scope], lengths={'chr': 9}, primary_CDS_chain_predictions=1)
        old = {**row, 'score': [-1.], 'score_kind': 'A6_complete_chain_logit_minus_window_null'}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); current = root/'new'; previous = root/'old'; current.mkdir(); previous.mkdir()
            (current/'scope.jsonl').write_text(json.dumps(row)+'\n')
            (previous/'scope.jsonl').write_text(json.dumps(old)+'\n')
            tx = as_transcript([(0, 9)], 'chr', '+', 'M28_A21_0_+_0', 1.)
            write_gff3(current/'A21.gff3', [tx], {'chr': 9}, 'M28_A21')
            with patch.object(V, 'INFER', current), patch.object(V.I.C, 'INFER', previous),\
                 patch.object(V.A, 'reference_context') as reference:
                predictions, _, _ = V.reconstruct(decoded, dict(scopes=[scope]))
                self.assertEqual(predictions['synthetic'], [tx]); reference.assert_not_called()
                with self.assertRaisesRegex(ValueError, 'prediction count'):
                    V.reconstruct({**decoded, 'primary_CDS_chain_predictions': 2}, dict(scopes=[scope]))
                reference.assert_not_called()

    def test_full_reference_and_background_denominators_are_not_reachable_subset(self):
        counts = V.H.A3.REF_COUNTS
        current = dict(per_species={s: dict(exact_chain=dict(reference=n), background=dict(bp=100))
                                   for s, n in counts.items()}, pooled=dict(exact_chain=dict(reference=7728)))
        V.H.check_ruler(current, current)
        bad = copy.deepcopy(current); species = next(iter(counts))
        bad['per_species'][species]['exact_chain']['reference'] -= 1
        with self.assertRaisesRegex(ValueError, 'denominators'):
            V.H.check_ruler(bad, current)
        bad = copy.deepcopy(current); bad['per_species'][species]['background']['bp'] = 99
        with self.assertRaisesRegex(ValueError, 'denominators'):
            V.H.check_ruler(bad, current)


if __name__ == '__main__':
    unittest.main()
