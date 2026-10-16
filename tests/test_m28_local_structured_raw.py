"""A44 paired raw interface on artificial inputs only."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import torch

import infer_local_structured_emission as I
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.export import write_gff3
from test_m28_emission_lattice import prefer_chain
from test_m28_local_structured_readout import fit_header


def receipts():
    fitted = {a: fit_header(a) for a in I.ARMS}; counts, refs = {}, {}
    for s, correct, nonmatch, unknown, long, eligible, five in (
            ('arabidopsis_thaliana', 1531, 347, 315, 548, 2172, 854),
            ('oryza_sativa', 478, 251, 132, 152, 744, 232)):
        old = dict(selected=correct+nonmatch+unknown, selected_correct=correct,
            selected_known_nonmatch=nonmatch, selected_unknown=unknown, selected_five_or_more_correct=long)
        counts[s] = dict(A31_A3=old)
        for delta, arm in enumerate(I.ARMS, 1):
            counts[s][arm] = dict(old, selected=old['selected']+delta,
                selected_correct=correct+delta, selected_five_or_more_correct=long+delta)
        refs[s] = dict(eligible=eligible, five_or_more=five)
    gates, supported = I.D.gates(counts)
    common = dict(completed=True, configuration=deepcopy(I.D.CONFIG), protocol=I.D.PROTOCOL,
        windows_per_arm=1491, original_draws=1536, original_candidates=16321, reference_injected_predictions=0,
        optimizer_steps_this_run=0, new_geometry_labels_computed=True, selected_counts_per_species=counts,
        eligible_reference_pairs=refs, predeclared_TRAIN_gates=gates, TRAIN_preparation_supported=supported,
        **{k: False for k in ('model_or_backbone_forward', 'new_candidate_search', 'new_backbone_forward',
            'R4_policy_or_catalog_changed', 'DEV_used', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened')})
    source = dict(fits=I.D.FIT_JOBS, potentials=I.POTENTIAL_JOB, decode=I.DECODE_JOB)
    prepared = dict(deepcopy(common), source_jobs=source,
        all_both_arm_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True)
    checked = dict(deepcopy(common), source_jobs=dict(source, readout=I.READOUT_JOB), candidate_window_pairs=19189,
                   **{k: True for k in I.CPU_PROOFS})
    return fitted, prepared, checked


class Emission:
    def __init__(self, offset=0.): self.offset, self.inputs = offset, []

    @property
    def query(self): raise AssertionError('No query generator allowed')

    def __call__(self, features, dna, valid):
        self.inputs.append((features, dna, valid))
        if tuple(dna.shape) != (1, I.M.WIDTH, 5) or valid != [features.shape[1]]:
            raise AssertionError('Original raw DNA width or valid length differs')
        return dict(segmentation_logits=features+self.offset)


class LoadedModel(torch.nn.Module):
    def __init__(self):
        super().__init__(); self.weight = torch.nn.Parameter(torch.zeros(366191))


class PairedRawTests(unittest.TestCase):
    def test_actual_two_arm_sources_seven_proofs_counts_and_both_comparisons_required(self):
        fitted, prepared, checked = receipts(); I.validate_preparation(fitted, prepared, checked)
        for proof in I.CPU_PROOFS:
            with self.subTest(proof=proof), self.assertRaisesRegex(ValueError, 'seven CPU'):
                I.validate_preparation(fitted, prepared, dict(checked, **{proof: False}))
        for change in ('fit_arm', 'source', 'baseline', 'denominator', 'scope', 'LOCAL_only', 'gates'):
            f, a, b = deepcopy(fitted), deepcopy(prepared), deepcopy(checked)
            if change == 'fit_arm': f['A44_LOCAL']['arm'] = 'A44_JOINT'
            elif change == 'source': b['source_jobs']['readout'] = '13424410'
            elif change == 'scope': a['DEV_used'] = True
            elif change == 'gates': a['predeclared_TRAIN_gates']['joint_vs_LOCAL']['pooled_five_or_more_correct_strictly_increased'] = False
            else:
                for report in (a, b):
                    if change == 'baseline': report['selected_counts_per_species']['arabidopsis_thaliana']['A31_A3']['selected_correct'] -= 1
                    elif change == 'denominator':
                        report['eligible_reference_pairs']['arabidopsis_thaliana']['eligible'] -= 1
                        report['eligible_reference_pairs']['oryza_sativa']['eligible'] += 1
                    else:
                        for values in report['selected_counts_per_species'].values(): values['A44_JOINT'] = deepcopy(values['A44_LOCAL'])
            with self.subTest(change=change), self.assertRaises(ValueError): I.validate_preparation(f, a, b)

    def test_failed_both_gate_missing_contract_or_explicit_input_precedes_models_and_FASTA(self):
        fitted, prepared, checked = receipts()
        with tempfile.TemporaryDirectory() as directory, patch.object(I, 'ROOT', Path(directory)),\
                patch.object(I, 'load_models') as models, patch.object(I.NATIVE, 'load_native') as native,\
                patch.object(I.NATIVE, 'input_contigs') as contigs:
            out = Path(directory)/'not_created'
            for failed in (True, False):
                def ready():
                    if failed:
                        a, b = deepcopy(prepared), deepcopy(checked)
                        for report in (a, b):
                            for values in report['selected_counts_per_species'].values():
                                values['A44_JOINT'] = deepcopy(values['A44_LOCAL'])
                            report['predeclared_TRAIN_gates'] = I.D.gates(report['selected_counts_per_species'])[0]
                        I.validate_preparation(fitted, a, b)
                    return fitted
                with patch.object(I, 'readiness', side_effect=ready), self.assertRaisesRegex(
                        ValueError, 'both fixed' if failed else 'separate fixed'):
                    I.main(out, Path(directory)/'unread.fa')
                self.assertFalse(out.exists())
            with patch.object(I, 'readiness', return_value=fitted), patch.object(Path, 'is_file', return_value=True),\
                    self.assertRaisesRegex(ValueError, 'Explicit authorized'):
                I.main(out, None)
            models.assert_not_called(); native.assert_not_called(); contigs.assert_not_called()

    def test_two_distinct_last_states_and_no_old_generator_load(self):
        models = [LoadedModel(), LoadedModel()]; fitted = receipts()[0]
        states = [dict(arm=a, configuration=I.D.F.CONFIG, protocol=I.D.F.PROTOCOL, step=4608,
            revision=I.D.R.C.I.F.REVISION, model={k: v.detach().clone() for k, v in m.state_dict().items()})
            for a, m in zip(I.ARMS, models)]
        with patch.object(torch, 'load', side_effect=states) as read,\
                patch.object(I.D.F.A.A36, 'new_model', side_effect=models) as create,\
                patch.object(I.D.R.C, 'load_generator', side_effect=AssertionError('Old weights forbidden')):
            loaded = I.load_models(fitted, 'cpu')
            self.assertEqual(create.call_count, 2); self.assertIsNot(loaded[I.ARMS[0]], loaded[I.ARMS[1]])
            self.assertEqual([c.args[0] for c in read.call_args_list], [I.ROOT/fitted[a]['final_checkpoint'] for a in I.ARMS])
            self.assertTrue(all(not m.training and all(not p.requires_grad for p in m.parameters()) for m in loaded.values()))
        with patch.object(torch, 'load', return_value=dict(states[0], arm='A44_JOINT')),\
                patch.object(I.D.F.A.A36, 'new_model') as create, self.assertRaises(ValueError):
            I.load_models(fitted, 'cpu')
        create.assert_not_called()

    def test_paired_stitch_exactly_matches_original_nonconstant_both_strand_context(self):
        sequence = ('ATGCCGTAACTG'*4)[:40]; locations = [[0, 20], [10, 30], [20, 40]]
        for strand in ('+', '-'):
            models = {a: Emission(float(10000*i)) for i, a in enumerate(I.ARMS)}
            tensors = [torch.arange(300, dtype=torch.float32).reshape(20, 15)+1000*i for i in range(3)]
            with patch.object(I.NATIVE, 'encode_window', side_effect=[v[None] for v in tensors]) as encoder:
                joined = I.paired_contig_emissions(models, None, None, sequence, locations, strand, 'cpu')
            self.assertEqual(encoder.call_count, 3)
            for left, right in zip(models[I.ARMS[0]].inputs, models[I.ARMS[1]].inputs):
                self.assertIs(left[0], right[0]); self.assertIs(left[1], right[1]); self.assertEqual(left[2], right[2])
            for j, arm in enumerate(I.ARMS):
                records = [dict(window_index=i, sequence=sequence[a:b] if strand == '+' else I.M.rc(sequence[a:b]),
                    local_potentials=tensors[i]+10000*j) for i, (a, b) in enumerate(locations)]
                original = stitch_contig_emissions(sequence, locations, strand, records)
                self.assertEqual(joined['sequence'], original['sequence'])
                self.assertEqual(joined['genomic_context_segments'], original['genomic_context_segments'])
                self.assertTrue(torch.equal(joined['potentials'][arm], original['local_potentials']))

    def test_normalized_minus_chain_phase_GFF_and_separate_arm_identity(self):
        oriented = 'ATGCGTTAGTTTAA'; raw = I.M.rc(oriented).lower()+'r'
        normalized, chain = 'N'+oriented, ((1, 5), (10, 15)); seen = []
        def encode(backbone, tokenizer, sequence, device):
            seen.append(sequence)
            return (prefer_chain(normalized, chain).float() if sequence == normalized else torch.zeros(15, 15))[None]
        models = {a: Emission() for a in I.ARMS}
        with tempfile.TemporaryDirectory() as directory, patch.object(I.NATIVE, 'encode_window', side_effect=encode):
            out = Path(directory)
            for arm in I.ARMS: (out/arm).mkdir()
            txs, scopes = I.annotate_contig(models, None, None, 'synthetic', 'chr', raw, 0, out, 'cpu')
            self.assertEqual(seen, [raw[:-1].upper()+'N', normalized])
            for arm in I.ARMS:
                self.assertEqual(len(models[arm].inputs), 2)
                minus = [t for t in txs[arm] if t['strand'] == '-']
                self.assertEqual(len(minus), 1); self.assertEqual(minus[0]['CDS'], [(0, 5, '2'), (10, 14, '0')])
                row = json.loads((out/scopes[arm][1]['file']).read_text())
                self.assertEqual(row['variant'], arm); self.assertEqual(row['source_fit_job'], I.D.FIT_JOBS[arm])
                self.assertEqual(row['checkpoint'], I.D.CONFIG['source_checkpoints'][arm]); self.assertEqual(row['chain_gains'], [140.])
                gff = out/arm/(arm+'.gff3'); write_gff3(gff, txs[arm], {'chr': 15}, 'M28_'+arm)
                cds = [l.split('\t') for l in gff.read_text().splitlines() if not l.startswith('#') and l.split('\t')[2] == 'CDS']
                self.assertEqual([(int(r[3]), int(r[4]), r[6], r[7]) for r in cds], [(1, 5, '-', '2'), (11, 14, '-', '0')])
                self.assertTrue(all(arm in t['id'] for t in txs[arm]))

    def test_empty_duplicate_streams_no_success_and_zero_prediction_receipt_keeps_both_arms(self):
        for rows in ([], [('user_input', 'dup', 'CCCCCC')]*2, [('user_input', 'one', 'CCCCCC')]):
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as directory,\
                    patch.object(I, 'ROOT', Path(directory)), patch.object(I, 'readiness', return_value=receipts()[0]),\
                    patch.object(Path, 'is_file', return_value=True),\
                    patch.object(I, 'load_models', return_value={a: LoadedModel().eval().requires_grad_(False) for a in I.ARMS}),\
                    patch.object(I.NATIVE, 'load_native', return_value=(torch.nn.Identity(), None)),\
                    patch.object(I.NATIVE, 'input_contigs', return_value=iter(rows)),\
                    patch.object(I, 'annotate_contig', return_value=({a: [] for a in I.ARMS}, {a: [] for a in I.ARMS})),\
                    patch.object(I, 'write_gff3') as gff, patch.object(I.D.R.C.R, 'write_summary') as summary:
                out = Path(directory)/'artificial'
                if len(rows) == 1:
                    I.main(out, Path(directory)/'not_opened.fa')
                    self.assertEqual(gff.call_count, 2); value = summary.call_args.args[1]
                    self.assertEqual(value['arms'], list(I.ARMS)); self.assertEqual(value['source_TRAIN_jobs'], I.source_jobs())
                    self.assertEqual(value['primary_CDS_chain_predictions'], {a: 0 for a in I.ARMS})
                    self.assertFalse(value['models_ensembled']); self.assertTrue(value['one_native_encoding_and_same_GEN_DNA_objects_per_paired_window'])
                else:
                    with self.assertRaisesRegex(ValueError, 'at least one' if not rows else 'unique'):
                        I.main(out, Path(directory)/'not_opened.fa')
                    gff.assert_not_called(); summary.assert_not_called()


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
