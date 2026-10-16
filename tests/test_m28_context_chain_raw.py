"""Synthetic A37 raw interface, frozen scientific gates and genomic coordinates."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_context_chain as I
import verify_context_chain_train as V
import evaluate_chain_quality as E


def receipts():
    fitted = dict(arm='A37', configuration=I.D.CONFIG, protocol=I.D.F.PROTOCOL,
        optimizer_steps=4608, feature_revision=I.D.F.A.R.C.I.F.REVISION, completed=True,
        preflight_job=I.D.F.PREFLIGHT_JOB, trainable_parameters=612324,
        initialization=I.D.CONFIG['initialization'], reference_injected_training_candidates=0,
        final_checkpoint=str((I.D.F.OUT/'step_004608.pt').relative_to(I.D.ROOT)),
        **{k:False for k in ('supervised_checkpoint_read','generator_forward','generator_updated',
            'new_candidate_search','new_backbone_forward','new_R4_labelling','DEV_used_for_selection',
            'test_setaria_used','A33_outputs_opened')})
    counts, refs = {}, {}
    for s, correct, wrong, five, eligible, long in (
            ('arabidopsis_thaliana',1531,347,548,2172,854), ('oryza_sativa',478,251,152,744,232)):
        old = dict(selected_correct=correct, selected_known_nonmatch=wrong, selected_five_or_more_correct=five)
        counts[s] = dict(A31_A3=old, A37={**old, 'selected_correct':correct+1, 'selected_five_or_more_correct':five+1})
        refs[s] = dict(eligible=eligible, five_or_more=long)
    common = dict(completed=True, configuration=I.D.CONFIG, protocol=I.D.PROTOCOL, windows=1491,
        natural_candidates=16321, reference_injected_predictions=0, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=I.D.gate(counts), TRAIN_preparation_supported=True)
    prepared = dict(**copy.deepcopy(common), source_fit_job=I.D.FIT_JOB, source_checks_job=I.CHECKS_JOB,
        checkpoint=fitted['final_checkpoint'], natural_score_file_closed_before_saved_label_arrays_open=True,
        original_natural_geometry_order_sources_unchanged=True, original_saved_labels_reused=True,
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        label_catalog_opened=False, new_backbone_forward=False)
    checked = dict(**copy.deepcopy(common), source_jobs=dict(fit=I.D.FIT_JOB, readout=I.READOUT_JOB),
        model_or_backbone_forward=False, independent_neural_forward_replay=False, **{k:True for k in I.CPU_PROOFS})
    return fitted, prepared, checked


class Generator:
    def __init__(self, query):
        self.query, self.inputs = query, []

    def __call__(self, features, dna, valid):
        self.inputs.append((features, dna, valid))
        return dict(features=features)

    def segmentation(self, hidden):
        raise AssertionError('Old A3 local emissions must not score A37 predictions')


class Scorer:
    def __init__(self, function):
        self.function, self.inputs = function, []

    def __call__(self, features, dna, valid, chains):
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw DNA, padding and valid-length contract changed')
        self.inputs.append((features, dna, valid, chains))
        return self.function(chains)


def manual_scores(chains, gain=None):
    exon = torch.full((len(chains),), -3.)
    junction = torch.tensor([float(len(c) > 1) for c in chains], dtype=torch.float32)
    additive = exon+junction; null = torch.tensor(2.)
    non = torch.tensor([32. if tuple(map(tuple,c)) == ((1,5),(10,15)) else 0. for c in chains], dtype=torch.float32)
    if gain is not None:
        non = torch.full((len(chains),), float(gain))+null-additive
    return dict(logits=additive+non-null, additive=additive, nonadditive=non,
                exon_additive=exon, junction_additive=junction, null=null,
                additive_only_logits=additive-null)


class ContextChainRawTests(unittest.TestCase):
    def test_actual_gate_cpu_identity_and_long_denominators_cannot_be_forged(self):
        fitted, prepared, checked = receipts(); I.validate_preparation(fitted, prepared, checked)
        for field, value in (('selected_correct',1531), ('selected_known_nonmatch',348)):
            bad, paired = copy.deepcopy(prepared), copy.deepcopy(checked)
            for report in (bad, paired):
                report['selected_counts_per_species']['arabidopsis_thaliana']['A37'][field] = value
            with self.assertRaisesRegex(ValueError, 'five-condition'):
                I.validate_preparation(fitted, bad, paired)
        wrong = {**checked, 'source_jobs':dict(fit=I.D.FIT_JOB, readout='wrong_job')}
        with self.assertRaisesRegex(ValueError, 'CPU reconciliation'):
            I.validate_preparation(fitted, prepared, wrong)
        bad, paired = copy.deepcopy(prepared), copy.deepcopy(checked)
        for report in (bad, paired):
            report['eligible_reference_pairs']['oryza_sativa']['five_or_more'] -= 1
        with self.assertRaisesRegex(ValueError, 'denominators'):
            I.validate_preparation(fitted, bad, paired)

    def test_failed_gate_or_missing_contract_prevent_native_model_fasta_and_output(self):
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

    def test_new_nonadditive_gain_shared_encoding_and_minus_gff_phase(self):
        oriented = 'ATGCGTTAGTTTAA'; query, features = fixture(oriented, [(0,4),(9,14)])
        generator, scorer = Generator(query), Scorer(manual_scores)
        raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            return (torch.cat([torch.zeros(1,4),features]) if sequence == 'N'+oriented else torch.zeros(15,4))[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', side_effect=encode):
            out = Path(directory)
            txs, scopes = I.annotate_contig(generator, scorer, None, None, 'synthetic', 'chr', raw, 0, out, 'cpu')
            self.assertEqual(seen, [raw[:-1].upper()+'N', 'N'+oriented])
            self.assertEqual(len(generator.inputs), 2); self.assertEqual(len(scorer.inputs), 2)
            for gen, new in zip(generator.inputs, scorer.inputs):
                self.assertIs(gen[0], new[0]); self.assertIs(gen[1], new[1]); self.assertEqual(gen[2], new[2])
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus), 1); self.assertEqual(minus[0]['CDS'], [(0,5,'2'),(10,14,'0')])
            row = next(I.C.I.saved_records(out/scopes[1]['file'])); V.verify_components({**row,'scores':row['score']})
            index = row['chains'].index([[1,5],[10,15]])
            self.assertEqual(row['score'][index], 28.); self.assertLess(row['additive'][index]-row['null'], 0.)
            self.assertEqual(minus[0]['score'], 28.); self.assertFalse(row['reference_used'])
            gff = out/'A37.gff3'; I.write_gff3(gff, txs, {'chr':15}, 'M28_A37')
            parsed = E.A.E.S.E.primary_transcripts(E.A.E.S.E.parse_annotation(gff, {'chr':15}), complete_only=False)
            E.A.validate_gff_roundtrip(parsed, txs)

    def test_nonpositive_empty_and_wrong_precision_keep_original_natural_contract(self):
        sequence = 'ATGAAATAA'; query, hidden = fixture(sequence, [(0,9)])
        generator = Generator(query); features = hidden[None]; dna = I.one_hot(sequence,I.M.WIDTH)[None]
        chains, sources = I.C.free_pool(generator, hidden, sequence)
        for gain in (0., -1.):
            scorer = Scorer(lambda c:manual_scores(c,gain))
            row = I.local_record(generator, hidden, scorer, features, dna, sequence, 0, 'synthetic')
            self.assertEqual(row['chains'], chains); self.assertEqual(row['decoder_sources'], sources)
            self.assertTrue(chains); self.assertTrue(all(v == gain for v in row['score']))
            V.verify_components({**row,'scores':row['score']})
            with tempfile.TemporaryDirectory() as directory, patch.object(I.C.I, 'encode_window', return_value=features):
                txs, _ = I.annotate_contig(generator, scorer, None, None, 'synthetic','chr',sequence,0,Path(directory),'cpu')
                self.assertEqual(txs, [])
        empty = I.local_record(generator, hidden, Scorer(manual_scores), features, dna, 'C'*len(sequence), 0, 'synthetic')
        self.assertEqual(empty['chains'], []); self.assertEqual(empty['score'], [])
        wrong = Scorer(lambda c:{k:v.double() for k,v in manual_scores(c).items()})
        with self.assertRaisesRegex(ValueError, 'FP32'):
            I.local_record(generator, hidden, wrong, features, dna, sequence, 0, 'synthetic')

    def test_empty_or_duplicate_explicit_contigs_cannot_write_complete_annotation(self):
        models = tuple(torch.nn.Linear(1,1).eval().requires_grad_(False) for _ in range(2))
        with tempfile.TemporaryDirectory() as directory, patch.object(I, 'ROOT', Path(directory)),\
             patch.object(I, 'readiness', return_value={'final_checkpoint':'unused'}),\
             patch.object(Path, 'is_file', return_value=True), patch.object(I, 'load_models', return_value=models),\
             patch.object(I.C.I, 'load_native', return_value=(None,None)),\
             patch.object(I.C.I, 'input_contigs') as contigs, patch.object(torch.cuda, 'reset_peak_memory_stats'),\
             patch.object(I, 'annotate_contig', return_value=([],[])) as annotate,\
             patch.object(I, 'write_gff3') as gff, patch.object(I.C.R, 'write_summary') as summary:
            for name, rows, error in (('empty',[], 'at least one'),
                    ('duplicate',[('user_input','chr','ATGAAATAA')]*2, 'unique')):
                contigs.return_value = iter(rows)
                with self.assertRaisesRegex(ValueError, error):
                    I.main(Path(directory)/name, Path('synthetic.fa'))
                self.assertFalse((Path(directory)/name/'summary.json').exists())
            self.assertEqual(annotate.call_count, 1); gff.assert_not_called(); summary.assert_not_called()


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
