"""Artificial A43 raw routing, joint context and export; no research inputs."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import torch

import infer_natural_competitor as I
from test_m28_emission_lattice import prefer_chain


def receipts():
    fitted = dict(arm='A43', configuration=deepcopy(I.D.F.CONFIG), protocol=I.D.F.PROTOCOL,
        optimizer_steps=4608, feature_revision=I.D.R.C.I.F.REVISION, completed=True,
        source_support_job='13424089', preflight_job='13424115', final_checkpoint=I.D.CONFIG['source_checkpoint'],
        trainable_parameters=366191, initialization='fresh_seed0_no_supervised_weight_inheritance',
        **{k: False for k in ('supervised_checkpoint_read', 'new_R4_labelling', 'generator_updated',
            'new_backbone_forward', 'DEV_used_for_selection', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened')})
    counts, refs = {}, {}
    for species, correct, nonmatch, unknown, long, eligible, five in (
            ('arabidopsis_thaliana', 1531, 347, 315, 548, 2172, 854),
            ('oryza_sativa', 478, 251, 132, 152, 744, 232)):
        old = dict(selected=correct+nonmatch+unknown, selected_correct=correct,
            selected_known_nonmatch=nonmatch, selected_unknown=unknown, selected_five_or_more_correct=long)
        counts[species] = dict(A31_A3=old, A43=dict(old, selected=old['selected']+1,
            selected_correct=correct+1, selected_five_or_more_correct=long+1))
        refs[species] = dict(eligible=eligible, five_or_more=five)
    common = dict(completed=True, configuration=deepcopy(I.D.CONFIG), protocol=I.D.PROTOCOL,
        windows=1491, original_draws=1536, original_candidates=16321, reference_injected_predictions=0,
        optimizer_steps_this_run=0, new_geometry_labels_computed=True,
        selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=I.D.gate(counts), TRAIN_preparation_supported=True,
        **{k: False for k in ('model_or_backbone_forward', 'new_candidate_search', 'R4_policy_or_catalog_changed',
            'DEV_used', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened')})
    source = dict(fit=I.D.FIT_JOB, potentials=I.POTENTIAL_JOB, decode=I.DECODE_JOB)
    prepared = dict(deepcopy(common), source_jobs=source,
        all_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True)
    checked = dict(deepcopy(common), source_jobs=dict(source, readout=I.READOUT_JOB), candidate_window_pairs=19189,
                   **{k: True for k in I.CPU_PROOFS})
    return fitted, prepared, checked


class Emission:
    """Artificial logits returned through the actual new single-model call."""
    def __init__(self): self.inputs = []

    @property
    def query(self): raise AssertionError('No query candidate generator allowed')

    def __call__(self, features, dna, valid):
        self.inputs.append((features, dna, valid))
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw DNA width or valid length differs')
        return dict(segmentation_logits=features)


class LoadedModel(torch.nn.Module):
    def __init__(self):
        super().__init__(); self.weight = torch.nn.Parameter(torch.zeros(366191))


class NaturalCompetitorRawTests(unittest.TestCase):
    def test_actual_A43_source_seven_proofs_baseline_and_full_denominators_required(self):
        fitted, prepared, checked = receipts(); I.validate_preparation(fitted, prepared, checked)
        for proof in I.CPU_PROOFS:
            with self.subTest(proof=proof), self.assertRaisesRegex(ValueError, 'seven CPU'):
                I.validate_preparation(fitted, prepared, dict(checked, **{proof: False}))
        for change in ('source', 'baseline', 'denominator', 'scope'):
            a, b = deepcopy(prepared), deepcopy(checked)
            if change == 'source': b['source_jobs']['readout'] = 'wrong'
            elif change == 'scope': a['DEV_used'] = True
            else:
                for report in (a, b):
                    if change == 'baseline': report['selected_counts_per_species']['arabidopsis_thaliana']['A31_A3']['selected_correct'] -= 1
                    else:
                        report['eligible_reference_pairs']['arabidopsis_thaliana']['eligible'] -= 1
                        report['eligible_reference_pairs']['oryza_sativa']['eligible'] += 1
            with self.subTest(change=change), self.assertRaises(ValueError): I.validate_preparation(fitted, a, b)

    def test_five_conditions_and_missing_contract_stop_before_native_models_or_FASTA(self):
        fitted, prepared, checked = receipts()
        with tempfile.TemporaryDirectory() as directory, patch.object(I, 'ROOT', Path(directory)),\
             patch.object(I, 'load_model') as model, patch.object(I.NATIVE, 'load_native') as native,\
             patch.object(I.NATIVE, 'input_contigs') as contigs:
            for failed in (True, False):
                def ready():
                    if failed:
                        a, b = deepcopy(prepared), deepcopy(checked)
                        for report in (a, b):
                            report['selected_counts_per_species']['arabidopsis_thaliana']['A43']['selected_correct'] = 1531
                        I.validate_preparation(fitted, a, b)
                    return fitted
                out = Path(directory)/'unexecuted'
                with patch.object(I, 'readiness', side_effect=ready), self.assertRaisesRegex(
                        ValueError, 'five-condition' if failed else 'separate fixed'):
                    I.main(out, Path(directory)/'not_read.fa')
                self.assertFalse(out.exists())
            with patch.object(I, 'readiness', return_value=fitted), patch.object(Path, 'is_file', return_value=True),\
                 self.assertRaisesRegex(ValueError, 'Explicit authorized FASTA'):
                I.main(Path(directory)/'unexecuted', None)
            model.assert_not_called(); native.assert_not_called(); contigs.assert_not_called()

    def test_last_A43_state_loads_one_frozen_model_without_old_generator(self):
        fitted = receipts()[0]; model = LoadedModel()
        state = dict(arm='A43', configuration=I.D.F.CONFIG, protocol=I.D.F.PROTOCOL, step=4608,
            revision=I.D.R.C.I.F.REVISION, model={k: v.detach().clone() for k, v in model.state_dict().items()})
        with patch.object(torch, 'load', return_value=state) as read,\
             patch.object(I.D.F.A.A36, 'new_model', return_value=model) as create,\
             patch.object(I.D.R.C, 'load_generator', side_effect=AssertionError('Old weights forbidden')):
            loaded = I.load_model(fitted, 'cpu')
            create.assert_called_once_with('cpu')
            self.assertEqual(read.call_args.args[0], I.ROOT/fitted['final_checkpoint'])
            self.assertIs(loaded, model); self.assertFalse(model.training)
            self.assertTrue(all(not p.requires_grad and p.grad is None for p in model.parameters()))
        with patch.object(torch, 'load', return_value=dict(state, step=1536)),\
             patch.object(I.D.F.A.A36, 'new_model') as create, self.assertRaises(ValueError):
            I.load_model(fitted, 'cpu')
        create.assert_not_called()

    def test_normalized_both_strands_single_emission_route_minus_phase_and_GFF(self):
        oriented = 'ATGCGTTAGTTTAA'; raw = I.M.rc(oriented).lower()+'r'
        normalized, chain = 'N'+oriented, ((1, 5), (10, 15))
        seen, model = [], Emission()
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            return (prefer_chain(normalized, chain).float() if sequence == normalized else torch.zeros(15, 15))[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.NATIVE, 'encode_window', side_effect=encode),\
             patch.object(I.D.R.C, 'free_pool', side_effect=AssertionError('No old candidate pool')):
            out = Path(directory)
            txs, scopes = I.annotate_contig(model, None, None, 'synthetic', 'chr', raw, 0, out, 'cpu')
            self.assertEqual(seen, [raw[:-1].upper()+'N', normalized]); self.assertEqual(len(model.inputs), 2)
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus), 1); self.assertEqual(minus[0]['CDS'], [(0, 5, '2'), (10, 14, '0')])
            row = json.loads((out/scopes[1]['file']).read_text())
            self.assertEqual(row['chain_gains'], [140.]); self.assertEqual(row['variant'], 'A43')
            self.assertFalse(row['reference_used']); self.assertEqual(row['context_policy'], I.CONTEXT_POLICY)
            gff = out/'toy.gff3'; I.write_gff3(gff, txs, {'chr': 15}, 'M28_A43')
            cds = [line.split('\t') for line in gff.read_text().splitlines() if not line.startswith('#') and line.split('\t')[2] == 'CDS']
            self.assertEqual([(int(r[3]), int(r[4]), r[6], r[7]) for r in cds], [(1, 5, '-', '2'), (11, 14, '-', '0')])

    def test_long_complete_chain_crosses_original_windows_before_joint_selection(self):
        sequence = 'C'*3000+'ATG'+'CCC'*9998+'TAA'+'C'*3864
        self.assertEqual(len(sequence), 36864); chain = ((3000, 33000),)
        logits = prefer_chain(sequence, chain).float(); model, seen = Emission(), []
        windows = [[0, 24576], [12288, 36864]]
        lookup = {sequence[a:b]: logits[a:b] for a, b in windows}
        def encode(backbone, tokenizer, oriented, device):
            seen.append(oriented)
            return lookup.get(oriented, torch.zeros(len(oriented), 15))[None]
        with tempfile.TemporaryDirectory() as directory, patch.object(I.NATIVE, 'encode_window', side_effect=encode):
            out = Path(directory)
            txs, scopes = I.annotate_contig(model, None, None, 'synthetic', 'long', sequence, 0, out, 'cpu')
            self.assertEqual(len(seen), 4); self.assertEqual(len(model.inputs), 4)
            self.assertEqual(scopes[0]['windows'], windows)
            self.assertFalse(any(a <= 3000 and 33000 <= b for a, b in windows))
            plus = [t for t in txs if t['strand'] == '+']
            self.assertEqual(len(plus), 1); self.assertEqual(plus[0]['CDS'], [(3000, 33000, '0')])
            row = json.loads((out/scopes[0]['file']).read_text())
            self.assertEqual(row['chain_gains'], [300000.]); self.assertEqual(row['gain'], 300000.)

    def test_zero_negative_no_ORF_and_empty_contig_keep_null_or_fail_without_encoder(self):
        for sequence, intergenic in (('ATGAAATAA', 0.), ('ATGAAATAA', 1.), ('CCCCCCC', 0.)):
            def encode(backbone, tokenizer, oriented, device):
                local = torch.zeros(len(oriented), 15); local[:, 0] = intergenic; return local[None]
            with self.subTest(sequence=sequence, intergenic=intergenic), tempfile.TemporaryDirectory() as directory,\
                 patch.object(I.NATIVE, 'encode_window', side_effect=encode):
                txs, scopes = I.annotate_contig(Emission(), None, None, 'synthetic', 'null', sequence, 0, Path(directory), 'cpu')
                self.assertEqual(txs, [])
                self.assertTrue(all(json.loads((Path(directory)/s['file']).read_text())['gain'] == 0 for s in scopes))
        with tempfile.TemporaryDirectory() as directory, patch.object(I.NATIVE, 'encode_window') as encoder,\
             self.assertRaisesRegex(ValueError, 'Empty contig'):
            I.annotate_contig(Emission(), None, None, 'synthetic', 'empty', '', 0, Path(directory), 'cpu')
        encoder.assert_not_called()

    def test_empty_and_duplicate_artificial_FASTA_streams_never_write_success(self):
        fitted = receipts()[0]
        for rows in ([], [('synthetic', 'dup', 'ATGTAA')]*2):
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as directory,\
                 patch.object(I, 'ROOT', Path(directory)), patch.object(I, 'readiness', return_value=fitted),\
                 patch.object(Path, 'is_file', return_value=True),\
                 patch.object(I, 'load_model', return_value=LoadedModel()),\
                 patch.object(I.NATIVE, 'load_native', return_value=(torch.nn.Identity(), None)),\
                 patch.object(I.NATIVE, 'input_contigs', return_value=iter(rows)),\
                 patch.object(I, 'annotate_contig', return_value=([], [])),\
                 patch.object(I, 'write_gff3') as gff, patch.object(I.D.R.C.R, 'write_summary') as summary:
                with self.assertRaisesRegex(ValueError, 'at least one' if not rows else 'unique'):
                    I.main(Path(directory)/'partial', Path(directory)/'not_read.fa')
                gff.assert_not_called(); summary.assert_not_called()


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
