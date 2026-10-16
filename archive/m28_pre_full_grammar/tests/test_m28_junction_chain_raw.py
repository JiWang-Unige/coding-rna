import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_junction_chain as I
from src.m28.junction_chain import JunctionChainHead
from src.m28.export import write_gff3
import evaluate_chain_quality as E


class SyntheticGenerator:
    def __init__(self, model):
        self.query = model; self.calls = 0

    def __call__(self, features, dna, valid):
        self.calls += 1
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw-window DNA/valid-length contract changed')
        return dict(features=features)


def synthetic_head(exon_gain=0., junction_gain=0.):
    head = JunctionChainHead(4).eval()
    with torch.no_grad():
        for parameter in head.parameters():
            parameter.zero_()
        head.additive.bias.fill_(exon_gain)
        head.junction_additive.bias.fill_(junction_gain)
    return head


def receipts():
    fitted = dict(arm='A21', completed=True, configuration=I.J.CONFIG, protocol=I.J.PROTOCOL,
        feature_revision=I.C.I.F.REVISION, optimizer_steps=4608,
        final_checkpoint=str((I.J.FIT/'step_004608.pt').relative_to(I.J.ROOT)),
        quality_trainable_parameters=248068, generator_trainable_parameters=0,
        epochs=[dict(steps=1536)]*3, DEV_used_for_selection=False, test_setaria_used=False)
    selected = {s: dict(A6=dict(selected_correct=4, selected_known_nonmatch=2, selected_five_or_more_correct=1),
                        A15=dict(selected_correct=5, selected_known_nonmatch=3, selected_five_or_more_correct=2),
                        A21=dict(selected_correct=5, selected_known_nonmatch=2, selected_five_or_more_correct=2))
                for s in I.C.R.ALLOW}
    common = dict(completed=True, configuration=I.J.CONFIG, protocol=I.J.PROTOCOL,
        natural_selection_per_species=selected, windows=1491, candidate_window_pairs=16321,
        optimizer_steps_this_run=0, DEV_forward_used=False, test_setaria_used=False,
        TRAIN_preparation_supported=True)
    prepared = dict(**common, checkpoint=fitted['final_checkpoint'],
        predeclared_TRAIN_preparation_gate=I.J.gate(selected), threshold=0., threshold_scan=False,
        new_candidate_search=False, new_backbone_forward=False, reference_injected_candidates=0,
        natural_candidate_geometry_labels_and_origins_unchanged=True, same_fit_immutable_hidden_used=True)
    checked = dict(**copy.deepcopy(common), original_frozen_TRAIN_gate=I.J.gate(selected),
        source_fit_readout_job='synthetic_fit', new_head_or_backbone_inference=False,
        actual_4608_steps_draw_weights_and_original_natural_domain_verified=True,
        actual_last_head_shapes_finite_tensors_components_selection_profiles_and_gate_verified=True,
        cache_zero_update_equivalence_verified=True, original_unknown_semantics_unchanged=True)
    return fitted, prepared, checked


class JunctionRawTests(unittest.TestCase):
    def test_recompute_gate_rejects_false_pass_missing_species_and_wrong_last_head(self):
        fitted, prepared, checked = receipts()
        I.validate_train_preparation(fitted, prepared, checked)
        species = next(iter(I.C.R.ALLOW))
        failed = copy.deepcopy(prepared)
        failed['natural_selection_per_species'][species]['A21']['selected_known_nonmatch'] = 3
        check_failed = {**checked, 'natural_selection_per_species': copy.deepcopy(failed['natural_selection_per_species'])}
        with self.assertRaisesRegex(ValueError, 'five-condition'):
            I.validate_train_preparation(fitted, failed, check_failed)
        incomplete = copy.deepcopy(prepared); incomplete['natural_selection_per_species'].pop(species)
        with self.assertRaisesRegex(ValueError, 'Both species'):
            I.validate_train_preparation(fitted, incomplete, checked)
        for key, value in (('optimizer_steps', 1536), ('final_checkpoint', 'early.pt'),
                           ('quality_trainable_parameters', 165635)):
            with self.subTest(key=key), self.assertRaises((ValueError, AssertionError)):
                I.validate_train_preparation({**fitted, key: value}, prepared, checked)

    def test_failed_gate_and_missing_dev_contract_reject_before_any_native_read_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); folders = [root/key for key in ('fit', 'readout', 'check')]
            for path in folders:
                path.mkdir(); (path/'STATUS').write_text('COMPLETED\n')
                (path/'JOBID').write_text('synthetic_fit\n')
            with patch.object(I, 'ROOT', root), patch.object(I.J, 'ROOT', root),\
                 patch.object(I.J, 'FIT', folders[0]), patch.object(I.J, 'READOUT', folders[1]),\
                 patch.object(I.J, 'CHECK', folders[2]), patch.object(I.torch, 'load') as load,\
                 patch.object(I.C, 'load_generator') as generator, patch.object(I.C.I, 'load_native') as native,\
                 patch.object(I.C.I, 'input_contigs') as contigs:
                for failed in (True, False):
                    fitted, prepared, checked = receipts()
                    if failed:
                        prepared['TRAIN_preparation_supported'] = False
                    for path, receipt in zip(folders, (fitted, prepared, checked)):
                        (path/'summary.json').write_text(json.dumps(receipt))
                    with self.assertRaisesRegex(ValueError, 'five-condition' if failed else 'separate fixed'):
                        I.main(root/'infer')
                    self.assertFalse((root/'infer').exists())
                for mocked in (load, generator, native, contigs):
                    mocked.assert_not_called()

    def test_real_junction_head_raw_normalization_minus_coordinates_and_gff_phase(self):
        oriented = 'ATGCGTTAGTTTAA'; model, features = fixture(oriented, [(0, 4), (9, 14)])
        generator = SyntheticGenerator(model); raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            values = torch.cat([torch.zeros(1, 4), features]) if sequence == 'N'+oriented else torch.zeros(15, 4)
            return values[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', side_effect=encode):
            out = Path(directory)
            txs, scopes = I.annotate_contig(generator, synthetic_head(1., 3.), None, None,
                                           'synthetic', 'chr', raw, 0, out, 'cpu')
            self.assertEqual(generator.calls, 2)
            self.assertEqual(seen, [raw[:-1].upper()+'N', 'N'+oriented])
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus), 1)
            self.assertEqual(minus[0]['CDS'], [(0, 5, '2'), (10, 14, '0')])
            row = next(I.C.I.saved_records(out/scopes[1]['file']))
            multi = next(i for i, chain in enumerate(row['chains']) if len(chain) == 2)
            self.assertEqual(row['score'][multi], 5.)
            self.assertEqual(row['score_components']['junction_additive'][multi], 3.)
            self.assertFalse(row['reference_used'])
            gff = out/'A21.gff3'; write_gff3(gff, txs, {'chr': 15}, 'M28_A21')
            parsed = E.A.E.S.E.primary_transcripts(E.A.E.S.E.parse_annotation(gff, {'chr': 15}), complete_only=False)
            E.A.validate_gff_roundtrip(parsed, txs)
            self.assertEqual([s['windows'] for s in scopes], [[[0, 15]], [[0, 15]]])

    def test_zero_gain_keeps_exact_free_pool_but_selects_null_and_negative_is_not_cut_early(self):
        sequence = 'ATGAAATAA'; model, features = fixture(sequence, [(0, 9)])
        generator = SyntheticGenerator(model)
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', return_value=features[None]):
            out = Path(directory)
            txs, scopes = I.annotate_contig(generator, synthetic_head(), None, None,
                                           'synthetic', 'chr', sequence, 0, out, 'cpu')
            self.assertEqual(txs, [])
            row = next(I.C.I.saved_records(out/scopes[0]['file']))
            expected, sources = I.C.free_pool(generator, features, sequence)
            serialized = [[list(exon) for exon in chain] for chain in expected]
            self.assertEqual(row['chains'], serialized); self.assertEqual(row['decoder_sources'], sources)
            self.assertTrue(row['chains']); self.assertTrue(all(v == 0. for v in row['score']))
            self.assertEqual(row['score_kind'], I.SCORE_KIND)
            self.assertEqual(scopes[0]['assembly_counts']['selected_unique'], 0)
            negative = I.quality_record(generator, synthetic_head(-2.), features, sequence, 0, 'synthetic|chr|0|+')
            self.assertEqual(negative['chains'], expected)
            self.assertTrue(all(v < 0 for v in negative['score']))


if __name__ == '__main__':
    unittest.main()
