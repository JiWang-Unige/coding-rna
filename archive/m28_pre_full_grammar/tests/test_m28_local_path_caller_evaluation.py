"""A31 full ruler and different-generator contrast, synthetic files only."""
from collections import Counter
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_local_path_caller as V
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def decoded_receipt(generator):
    scopes = []
    for ordinal, (species, lengths) in enumerate(V.I.C.I.DEV_LENGTHS.items()):
        for seqid, length in lengths.items():
            for strand in ('+', '-'):
                scopes.append(dict(species=species, seqid=seqid, strand=strand, ordinal=ordinal,
                    windows=[[a, min(a+V.I.M.WIDTH, length)] for a in V.I.M.starts(length)],
                    file='contig_'+str(ordinal)+'_'+strand+'.jsonl'))
    return dict(arm='A31', variant='A31_'+generator, completed=True, inference_complete=True,
        configuration=V.I.R.CONFIG, protocol=V.I.DEV_PROTOCOL, TRAIN_protocol=V.I.R.PROTOCOL,
        TRAIN_preparation_gate_passed=True, generator_checkpoint=V.I.R.CONFIG['generators'][generator],
        generator_optimizer_steps=4608, new_parameters=0, feature_revision=V.I.C.I.F.REVISION,
        feature_spec=V.I.C.R.FEATURE_SPEC, input_mode='fixed_DEV_raw_FASTA', reference_used=False,
        labels_or_feature_cache_used=False, optimizer_steps_this_run=0, DEV_used_for_selection=False,
        test_setaria_used=False, score_kind=V.I.SCORE_KIND, windows=8690, scopes=scopes,
        lengths={q: n for qs in V.I.C.I.DEV_LENGTHS.values() for q, n in qs.items()})


def raw_row(generator='A3'):
    return dict(window_index=0, window_id='synthetic|chr|0|+', chains=[[[0, 9]]], score=[1.],
        decoder_sources=[['beam4', 'greedy']], score_kind=V.I.SCORE_KIND, variant='A31_'+generator,
        score_threshold=0., reference_used=False, counts=dict(free_candidates=1, nonpositive_candidates=0))


def old_row():
    return dict(**{k: v for k, v in raw_row().items() if k not in ('score', 'score_kind')},
        score=[-1.], score_kind='A6_complete_chain_logit_minus_window_null', quality_logit_threshold=0.)


class LocalPathEvaluationTests(unittest.TestCase):
    def test_full_grid_context_frozen_model_and_variant_required(self):
        for generator in ('A3', 'A23'):
            good = decoded_receipt(generator); V.validate_complete(good, generator)
            for key, value in (('scopes', good['scopes'][:-1]), ('input_mode', 'explicit_raw_FASTA'),
                               ('feature_spec', {'encoder_context': 'whole_4096'}), ('reference_used', True),
                               ('generator_checkpoint', 'early.pt'), ('generator_optimizer_steps', 1536)):
                with self.subTest(generator=generator, key=key), self.assertRaises(ValueError):
                    V.validate_complete({**good, key: value}, generator)
            truncated = copy.deepcopy(good); truncated['scopes'][0]['windows'].pop()
            with self.assertRaises(ValueError):
                V.validate_complete(truncated, generator)

    def test_A3_requires_same_pool_A23_allows_different_pool_but_same_input(self):
        scope = dict(species='synthetic', seqid='chr', strand='+', windows=[[0, 9]])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); current, prior = root/'new.jsonl', root/'old.jsonl'
            prior.write_text(json.dumps(old_row())+'\n'); current.write_text(json.dumps(raw_row())+'\n')
            self.assertEqual(len(list(V.paired_records(current, scope, prior, scope, 'A3', Counter()))), 1)
            changed = {**raw_row(), 'chains': [[[0, 6]]]}
            current.write_text(json.dumps(changed)+'\n')
            with self.assertRaisesRegex(ValueError, 'Same frozen A3'):
                list(V.paired_records(current, scope, prior, scope, 'A3', Counter()))
            changed['variant'] = 'A31_A23'; current.write_text(json.dumps(changed)+'\n')
            self.assertEqual(len(list(V.paired_records(current, scope, prior, scope, 'A23', Counter()))), 1)
            for key, value in (('window_id', 'synthetic|chr|1|+'), ('score', [float('nan')]),
                               ('decoder_sources', [['invented']]), ('score_threshold', .2),
                               ('counts', dict(free_candidates=1, nonpositive_candidates=1))):
                current.write_text(json.dumps({**changed, key: value})+'\n')
                with self.subTest(key=key), self.assertRaises(ValueError):
                    list(V.paired_records(current, scope, prior, scope, 'A23', Counter()))
            current.write_text('')
            with self.assertRaises(ValueError):
                list(V.paired_records(current, scope, prior, scope, 'A23', Counter()))

    def test_saved_selection_GFF_phase_and_count_close_before_reference(self):
        row = raw_row()
        scope = dict(species='synthetic', seqid='chr', strand='+', ordinal=0, windows=[[0, 9]], file='scope.jsonl',
            candidate_counts=row['counts'], assembly_counts=assemble_predictions([row], [[0, 9]], '+')['counts'])
        decoded = dict(scopes=[scope], lengths={'chr': 9}, primary_CDS_chain_predictions=1)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); current, prior = root/'new', root/'old'; current.mkdir(); prior.mkdir()
            (current/'scope.jsonl').write_text(json.dumps(row)+'\n')
            (prior/'scope.jsonl').write_text(json.dumps(old_row())+'\n')
            tx = as_transcript([(0, 9)], 'chr', '+', 'M28_A31_A3_0_+_0', 1.)
            write_gff3(current/'A31_A3.gff3', [tx], {'chr': 9}, 'M28_A31_A3')
            with patch.object(V.I.C, 'INFER', prior), patch.object(V.A, 'reference_context') as reference:
                predictions, _, _ = V.reconstruct(decoded, dict(scopes=[scope]), current, 'A3')
                self.assertEqual(predictions['synthetic'], [tx]); reference.assert_not_called()
                with self.assertRaisesRegex(ValueError, 'prediction count'):
                    V.reconstruct({**decoded, 'primary_CDS_chain_predictions': 2}, dict(scopes=[scope]), current, 'A3')
                reference.assert_not_called()

    def test_complete_reference_and_background_denominators_not_reachable_subset(self):
        current = dict(per_species={s: dict(exact_chain=dict(reference=n), background=dict(bp=100))
            for s, n in V.H.A3.REF_COUNTS.items()}, pooled=dict(exact_chain=dict(reference=7728)))
        V.H.check_ruler(current, current)
        bad = copy.deepcopy(current); species = next(iter(V.H.A3.REF_COUNTS))
        bad['per_species'][species]['exact_chain']['reference'] -= 1
        with self.assertRaisesRegex(ValueError, 'denominators'):
            V.H.check_ruler(bad, bad)
        bad = copy.deepcopy(current); bad['per_species'][species]['background']['bp'] -= 1
        with self.assertRaisesRegex(ValueError, 'denominators'):
            V.H.check_ruler(bad, current)


if __name__ == '__main__':
    unittest.main()
