"""Synthetic A31 raw interface, scientific preconditions and GFF roundtrip."""
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
import infer_local_path_caller as I
from src.m28.phase_evidence_chain import candidate_roles
from src.m28.export import write_gff3
import evaluate_chain_quality as E


class SyntheticGenerator:
    def __init__(self, query, segmentation):
        self.query, self.segmentation, self.calls = query, segmentation, 0

    def __call__(self, features, dna, valid):
        self.calls += 1
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw DNA and valid-length contract changed')
        return dict(features=features)


def receipts():
    counts = {s: {a: dict(selected_correct=6 if a.startswith('A31') else 4,
        selected_known_nonmatch=2, selected_five_or_more_correct=2 if a.startswith('A31') else 1)
        for a in I.R.ARMS} for s in I.C.R.ALLOW}
    references = {s: dict(eligible=10, five_or_more=3) for s in I.C.R.ALLOW}
    common = dict(arm='A31', completed=True, configuration=I.R.CONFIG, protocol=I.R.PROTOCOL,
        windows=1491, optimizer_steps_this_run=0, reference_injected_candidates=0,
        new_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False,
        test_setaria_used=False, selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=I.R.gate(counts, references),
        TRAIN_preparation_supported={'A31_A3': True, 'A31_A23': True})
    prepared = dict(**copy.deepcopy(common), natural_score_file_closed_before_catalog_open=True,
        original_greedy_source_replay_verified=True, model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True, original_natural_pools_unchanged=True)
    checked = dict(**copy.deepcopy(common), source_jobs={'readout': 'synthetic_readout'},
        actual_saved_natural_pool_and_greedy_subset_identity_verified=True,
        actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed=True,
        actual_five_fixed_gate_integer_arithmetic_verified=True,
        original_last_checkpoint_tensors_finite_and_parameter_counts_verified=True,
        new_head_or_backbone_forward=False)
    return prepared, checked


class LocalPathRawTests(unittest.TestCase):
    def test_independent_variant_gate_recomputed_and_missing_species_rejected(self):
        prepared, checked = receipts()
        for a in ('A3', 'A23'):
            I.validate_train_preparation(prepared, checked, a)
        species = next(iter(I.C.R.ALLOW))
        bad = copy.deepcopy(prepared); bad['selected_counts_per_species'][species]['A31_A3']['selected_known_nonmatch'] = 9
        paired = {**checked, 'selected_counts_per_species': copy.deepcopy(bad['selected_counts_per_species'])}
        with self.assertRaisesRegex(ValueError, 'five-condition'):
            I.validate_train_preparation(bad, paired, 'A3')
        # One failed variant does not close a different variant with reconciled gates.
        actual = I.R.gate(bad['selected_counts_per_species'], bad['eligible_reference_pairs'])
        for rec in (bad, paired):
            rec['predeclared_TRAIN_gate'] = actual
            rec['TRAIN_preparation_supported'] = {a: all(v.values()) for a, v in actual.items()}
        I.validate_train_preparation(bad, paired, 'A23')
        incomplete = copy.deepcopy(prepared); incomplete['eligible_reference_pairs'].pop(species)
        with self.assertRaisesRegex(ValueError, 'Both species'):
            I.validate_train_preparation(incomplete, checked, 'A3')

    def test_failed_gate_or_missing_contract_does_not_open_native_model_fasta_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); readout = root/'readout'; readout.mkdir()
            (readout/'JOBID').write_text('synthetic_readout\n')
            with patch.object(I, 'ROOT', root), patch.object(I.R, 'READOUT', readout),\
                 patch.object(I.R, 'completed') as completed, patch.object(I, 'load_generator') as generator,\
                 patch.object(I.C.I, 'load_native') as native, patch.object(I.C.I, 'input_contigs') as contigs:
                for failed in (True, False):
                    prepared, checked = receipts()
                    if failed:
                        prepared['TRAIN_preparation_supported']['A31_A3'] = False
                    completed.side_effect = [prepared, checked]
                    with self.assertRaisesRegex(ValueError, 'five-condition' if failed else 'separate fixed'):
                        I.main(root/'infer', 'A3')
                    self.assertFalse((root/'infer').exists())
                for mocked in (generator, native, contigs):
                    mocked.assert_not_called()

    def test_minus_coordinates_raw_normalization_phase_and_gff_roundtrip(self):
        oriented = 'ATGCGTTAGTTTAA'; query, features = fixture(oriented, [(0, 4), (9, 14)])
        def emission(hidden):
            logits = torch.zeros(len(hidden), 15)
            if bool(hidden.any()):
                logits.fill_(-2); logits[:, 0] = 0
                roles = candidate_roles(((1, 5), (10, 15)), len(hidden), hidden.device)
                logits[torch.arange(len(hidden)), roles] = 2
            return logits
        generator = SyntheticGenerator(query, emission)
        raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            return (torch.cat([torch.zeros(1, 4), features]) if sequence == 'N'+oriented else torch.zeros(15, 4))[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', side_effect=encode):
            out = Path(directory)
            txs, scopes = I.annotate_contig(generator, None, None, 'synthetic', 'chr', raw, 0, out, 'cpu', 'A3')
            self.assertEqual(seen, [raw[:-1].upper()+'N', 'N'+oriented])
            self.assertEqual(generator.calls, 2)
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus), 1)
            self.assertEqual(minus[0]['CDS'], [(0, 5, '2'), (10, 14, '0')])
            row = next(I.C.I.saved_records(out/scopes[1]['file']))
            multi = row['chains'].index([[1, 5], [10, 15]])
            self.assertEqual(row['score'][multi], 28.)
            self.assertFalse(row['reference_used'])
            gff = out/'A31.gff3'; write_gff3(gff, txs, {'chr': 15}, 'M28_A31_A3')
            parsed = E.A.E.S.E.primary_transcripts(E.A.E.S.E.parse_annotation(gff, {'chr': 15}), complete_only=False)
            E.A.validate_gff_roundtrip(parsed, txs)

    def test_zero_selects_null_negative_retains_natural_pool_and_empty_pool_is_supported(self):
        sequence = 'ATGAAATAA'; query, hidden = fixture(sequence, [(0, 9)])
        generator = SyntheticGenerator(query, lambda h: torch.zeros(len(h), 15))
        expected, sources = I.C.free_pool(generator, hidden, sequence)
        zero = I.local_record(generator, hidden, sequence, 0, 'synthetic|chr|0|+', 'A23')
        self.assertEqual(zero['chains'], expected); self.assertEqual(zero['decoder_sources'], sources)
        self.assertTrue(zero['chains']); self.assertTrue(all(v == 0 for v in zero['score']))
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', return_value=hidden[None]):
            txs, scopes = I.annotate_contig(generator, None, None, 'synthetic', 'chr', sequence, 0, Path(directory), 'cpu', 'A23')
            self.assertEqual(txs, [])
            self.assertTrue(all(s['assembly_counts']['selected_unique'] == 0 for s in scopes))
        generator.segmentation = lambda h: torch.cat([torch.zeros(len(h), 1), -torch.ones(len(h), 14)], -1)
        negative = I.local_record(generator, hidden, sequence, 0, 'synthetic|chr|0|+', 'A3')
        self.assertEqual(negative['chains'], expected); self.assertTrue(all(v < 0 for v in negative['score']))
        empty = I.local_record(generator, hidden, 'C'*len(sequence), 0, 'synthetic|chr|0|+', 'A3')
        self.assertEqual(empty['chains'], []); self.assertEqual(empty['score'], [])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
