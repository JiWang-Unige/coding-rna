"""Synthetic A36 two-representation caller, unchanged gates and GFF geometry."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_structured_emission as I
import evaluate_chain_quality as E
from src.m28.phase_evidence_chain import candidate_roles


def receipts():
    fitted = dict(arm='A36', configuration=I.D.CONFIG, protocol=I.D.F.PROTOCOL,
        optimizer_steps=4608, feature_revision=I.D.F.A.C.I.F.REVISION, completed=True,
        preflight_job=I.D.F.PREFLIGHT_JOB, final_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.D.ROOT)),
        trainable_parameters=366191, DEV_used_for_selection=False, test_setaria_used=False, A33_outputs_opened=False)
    counts = {s: dict(A31_A3=dict(selected_correct=10, selected_known_nonmatch=3, selected_five_or_more_correct=2),
        A36=dict(selected_correct=11, selected_known_nonmatch=3, selected_five_or_more_correct=3)) for s in I.C.R.ALLOW}
    common = dict(completed=True, configuration=I.D.CONFIG, protocol=I.D.PROTOCOL, windows=1491,
        natural_candidates=16321, reference_injected_predictions=0, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, selected_counts_per_species=counts,
        eligible_reference_pairs={s:dict(eligible=1458, five_or_more=3) for s in counts},
        predeclared_TRAIN_gate=I.D.gate(counts), TRAIN_preparation_supported=True)
    prepared = dict(**copy.deepcopy(common), source_fit_job=I.D.FIT_JOB, checkpoint=fitted['final_checkpoint'],
        natural_score_file_closed_before_saved_label_arrays_open=True, original_natural_geometry_order_sources_unchanged=True,
        original_saved_labels_reused=True, model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        label_catalog_opened=False, new_backbone_forward=False)
    checked = dict(**copy.deepcopy(common), source_jobs=dict(fit=I.D.FIT_JOB, readout=I.READOUT_JOB),
        model_or_backbone_forward=False, **{k:True for k in I.CPU_PROOFS})
    return fitted, prepared, checked


class Generator:
    def __init__(self, query):
        self.query, self.inputs = query, []

    def __call__(self, features, dna, valid):
        self.inputs.append((features, dna, valid))
        return dict(features=features)

    def segmentation(self, hidden):
        raise AssertionError('Old A3 emissions must not score A36 predictions')


class Emission:
    def __init__(self, function):
        self.function, self.inputs = function, []

    def __call__(self, features, dna, valid):
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw DNA, padding and valid-length contract changed')
        self.inputs.append((features, dna, valid))
        return dict(segmentation_logits=self.function(features[0])[None])


class StructuredEmissionRawTests(unittest.TestCase):
    def test_both_species_stronger_gate_and_cpu_identity_cannot_be_forged(self):
        fitted, prepared, checked = receipts(); I.validate_preparation(fitted, prepared, checked)
        s = next(iter(I.C.R.ALLOW))
        for field, value in (('selected_correct', 10), ('selected_known_nonmatch', 4)):
            bad, paired = copy.deepcopy(prepared), copy.deepcopy(checked)
            for report in (bad, paired):
                report['selected_counts_per_species'][s]['A36'][field] = value
            with self.assertRaisesRegex(ValueError, 'five-condition'):
                I.validate_preparation(fitted, bad, paired)
        wrong = {**checked, 'source_jobs':dict(fit=I.D.FIT_JOB, readout='wrong_job')}
        with self.assertRaisesRegex(ValueError, 'CPU reconciliation'):
            I.validate_preparation(fitted, prepared, wrong)

    def test_failed_gate_and_missing_contract_prevent_native_weights_fasta_and_output(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(I, 'ROOT', Path(directory)),\
             patch.object(I.D.F.A, 'require_jobs'), patch.object(I.D.R, 'completed') as completed,\
             patch.object(I, 'load_models') as models, patch.object(I.C.I, 'load_native') as native,\
             patch.object(I.C.I, 'input_contigs') as contigs:
            for failed in (True, False):
                fitted, prepared, checked = receipts()
                if failed:
                    prepared['TRAIN_preparation_supported'] = False
                completed.side_effect = [fitted, prepared, checked]
                with self.assertRaisesRegex(ValueError, 'five-condition' if failed else 'separate fixed'):
                    I.main(Path(directory)/'native')
                self.assertFalse((Path(directory)/'native').exists())
            for mocked in (models, native, contigs):
                mocked.assert_not_called()

    def test_new_emission_scores_share_encoding_and_preserve_minus_gff_phase(self):
        oriented = 'ATGCGTTAGTTTAA'; query, features = fixture(oriented, [(0, 4), (9, 14)])
        def emit(hidden):
            logits = torch.zeros(len(hidden), 15)
            if bool(hidden.any()):
                logits.fill_(-2); logits[:, 0] = 0
                roles = candidate_roles(((1, 5), (10, 15)), len(hidden))
                logits[torch.arange(len(hidden)), roles] = 2
            return logits
        generator, emission = Generator(query), Emission(emit)
        raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            return (torch.cat([torch.zeros(1, 4), features]) if sequence == 'N'+oriented else torch.zeros(15, 4))[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', side_effect=encode):
            out = Path(directory)
            txs, scopes = I.annotate_contig(generator, emission, None, None, 'synthetic', 'chr', raw, 0, out, 'cpu')
            self.assertEqual(seen, [raw[:-1].upper()+'N', 'N'+oriented])
            self.assertEqual(len(generator.inputs), 2); self.assertEqual(len(emission.inputs), 2)
            for gen, seg in zip(generator.inputs, emission.inputs):
                self.assertIs(gen[0], seg[0]); self.assertIs(gen[1], seg[1]); self.assertEqual(gen[2], seg[2])
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus), 1); self.assertEqual(minus[0]['CDS'], [(0, 5, '2'), (10, 14, '0')])
            row = next(I.C.I.saved_records(out/scopes[1]['file']))
            self.assertEqual(row['score'][row['chains'].index([[1, 5], [10, 15]])], 28.)
            self.assertFalse(row['reference_used'])
            gff = out/'A36.gff3'; I.write_gff3(gff, txs, {'chr':15}, 'M28_A36')
            parsed = E.A.E.S.E.primary_transcripts(E.A.E.S.E.parse_annotation(gff, {'chr':15}), complete_only=False)
            E.A.validate_gff_roundtrip(parsed, txs)

    def test_zero_negative_and_empty_natural_pools_do_not_change_proposal_geometry(self):
        sequence = 'ATGAAATAA'; query, hidden = fixture(sequence, [(0, 9)])
        generator = Generator(query); zero = torch.zeros(len(sequence), 15)
        chains, sources = I.C.free_pool(generator, hidden, sequence)
        record = I.local_record(generator, hidden, zero, sequence, 0, 'synthetic')
        self.assertEqual(record['chains'], chains); self.assertEqual(record['decoder_sources'], sources)
        self.assertTrue(chains); self.assertTrue(all(v == 0 for v in record['score']))
        emission = Emission(lambda h:torch.zeros(len(h), 15))
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', return_value=hidden[None]):
            txs, _ = I.annotate_contig(generator, emission, None, None, 'synthetic', 'chr', sequence, 0, Path(directory), 'cpu')
            self.assertEqual(txs, [])
        negative = zero.clone(); negative[:, 0] = 1
        row = I.local_record(generator, hidden, negative, sequence, 0, 'synthetic')
        self.assertEqual(row['chains'], chains); self.assertTrue(all(v < 0 for v in row['score']))
        empty = I.local_record(generator, hidden, zero, 'C'*len(sequence), 0, 'synthetic')
        self.assertEqual(empty['chains'], []); self.assertEqual(empty['score'], [])


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
