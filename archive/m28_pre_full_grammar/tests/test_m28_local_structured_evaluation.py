"""A44 paired raw evaluator on artificial JSON/GFF and references only."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import evaluate_local_structured_emission as W
from src.m28.export import as_transcript, write_gff3
from test_m28_natural_competitor_evaluation import full_context, comparators


def saved_pair(root, contigs, calls):
    scopes, totals = {}, {}
    for arm in W.I.ARMS:
        (root/arm).mkdir(); scopes[arm], txs = [], []
        for ordinal, (seqid, c) in enumerate(contigs.items()):
            length = c['length']; windows = [[a, min(a+W.I.M.WIDTH, length)] for a in W.I.M.starts(length)]
            for strand, tag in (('+', 'p'), ('-', 'm')):
                chains = [[list(p) for p in chain] for chain in calls.get((arm, seqid, strand), [])]
                gains = [float(i+1) for i in range(len(chains))]
                converted = [[list(p) for p in W.genomic_chain(c, 0, length, strand)] for c in chains]
                name = arm+'/contig_'+str(ordinal).zfill(4)+'_'+tag+'.json'
                scopes[arm].append(dict(species='user_input', seqid=seqid, strand=strand, ordinal=ordinal,
                    arm=arm, windows=windows, file=name, context_policy=W.I.CONTEXT_POLICY,
                    primary_CDS_chain_predictions=len(chains), decode_seconds=0., shared_paired_encoding_seconds=0.))
                row = dict(variant=arm, score_kind=W.I.SCORE_KIND, species='user_input', seqid=seqid,
                    strand=strand, length=length, windows=windows, source_fit_job=W.I.D.FIT_JOBS[arm],
                    checkpoint=W.I.D.CONFIG['source_checkpoints'][arm], chains=chains, chain_gains=gains,
                    gain=sum(gains), genomic_chains=converted, counts={}, reference_used=False,
                    reference_used_in_scoring=False, context_policy=W.I.CONTEXT_POLICY,
                    genomic_context_segments=[list(s) for s in W.context_segments(length, windows)])
                (root/name).write_text(json.dumps(row)+'\n')
                for i, (chain, gain) in enumerate(zip(converted, gains)):
                    txs.append(as_transcript(chain, seqid, strand, 'M28_'+arm+'_'+str(ordinal)+'_'+tag+'_'+str(i), gain))
        totals[arm] = len(txs)
        write_gff3(root/arm/(arm+'.gff3'), txs, {q: c['length'] for q, c in contigs.items()}, 'M28_'+arm)
    return dict(arms=list(W.I.ARMS), completed=True, inference_complete=True,
        TRAIN_configuration=W.I.D.CONFIG, protocol=W.I.RAW_PROTOCOL, TRAIN_protocol=W.I.D.PROTOCOL,
        source_TRAIN_jobs=W.I.source_jobs(), checkpoints=W.I.D.CONFIG['source_checkpoints'],
        optimizer_steps_source_per_arm=4608, model_parameters_per_arm=366191,
        feature_revision=W.I.NATIVE.F.REVISION, feature_spec=W.I.D.R.C.R.FEATURE_SPEC,
        input_mode='explicit_raw_FASTA', input_fasta='/synthetic/not_opened.fa',
        lengths={q: c['length'] for q, c in contigs.items()}, scopes=scopes,
        context_policy=W.I.CONTEXT_POLICY, score_kind=W.I.SCORE_KIND, max_exons=128, threshold=0., null=0.,
        optimizer_steps_this_run=0, new_parameters=0, primary_CDS_chain_predictions=totals,
        **{k: True for k in ('one_native_encoding_and_same_GEN_DNA_objects_per_paired_window',
            'frozen_source_weights_unchanged', 'all_source_parameters_frozen_and_gradients_absent')},
        **{k: False for k in ('reference_used', 'labels_or_feature_cache_used', 'generator_model_loaded',
            'query_or_beam_decoding', 'local_selected_only_owner_assembly', 'models_ensembled', 'pretrained_encoder_updated')})


class PairedEvaluationTests(unittest.TestCase):
    def test_two_arm_source_receipt_checkpoint_and_scope_identity(self):
        contigs = dict(toy=dict(species='synthetic', length=36864))
        with tempfile.TemporaryDirectory() as directory:
            decoded = saved_pair(Path(directory), contigs, {}); W.validate_complete(decoded, contigs)
            for key, wrong in (('arms', ['A44_JOINT']), ('model_parameters_per_arm', 649457),
                    ('input_mode', 'fixed_DEV_raw_FASTA'), ('models_ensembled', True),
                    ('source_TRAIN_jobs', {'fit': '13424143'}), ('threshold', .1),
                    ('one_native_encoding_and_same_GEN_DNA_objects_per_paired_window', False)):
                with self.subTest(key=key), self.assertRaises(ValueError): W.validate_complete(dict(decoded, **{key: wrong}), contigs)
            for change in ('missing', 'swapped', 'grid', 'context', 'checkpoint'):
                bad = deepcopy(decoded)
                if change == 'missing': bad['scopes']['A44_LOCAL'].pop()
                elif change == 'swapped': bad['scopes']['A44_JOINT'][0]['arm'] = 'A44_LOCAL'
                elif change == 'grid': bad['scopes']['A44_LOCAL'][0]['windows'][1][0] += 1
                elif change == 'context': bad['scopes']['A44_LOCAL'][0]['context_policy'] = 'old_chain_owner'
                else: bad['checkpoints']['A44_LOCAL'] = bad['checkpoints']['A44_JOINT']
                with self.subTest(change=change), self.assertRaises(ValueError): W.validate_complete(bad, contigs)

    def test_minus_phase_gain_pairing_and_cross_window_chain_in_separate_GFF(self):
        contigs = dict(toy=dict(species='synthetic', length=36864))
        calls = {('A44_LOCAL', 'toy', '-'): [[(0, 4), (9, 14)], [(20, 26)]],
                 ('A44_JOINT', 'toy', '+'): [[(3000, 33000)]]}
        with tempfile.TemporaryDirectory() as directory, patch.object(W.V.A, 'reference_context') as reference:
            root = Path(directory); decoded = saved_pair(root, contigs, calls)
            predictions, counts = W.reconstruct(decoded, root, contigs)
            self.assertEqual(predictions['A44_LOCAL']['synthetic'][0]['CDS'], [(36850, 36855, '2'), (36860, 36864, '0')])
            self.assertEqual([t['score'] for t in predictions['A44_LOCAL']['synthetic']], [1., 2.])
            self.assertEqual(predictions['A44_JOINT']['synthetic'][0]['CDS'], [(3000, 33000, '0')])
            self.assertEqual(counts['A44_LOCAL']['toy|+']['selected_chains'], 0)
            self.assertTrue(all(arm in t['id'] for arm in W.I.ARMS for t in predictions[arm]['synthetic']))
            reference.assert_not_called()

    def test_any_arm_bad_saved_source_gain_GFF_or_count_stops_before_reference(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        calls = {(arm, first, '+'): [[(0, 6)]] for arm in W.I.ARMS}
        for change in ('source', 'gain', 'GFF_phase', 'count'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory, patch.object(W.V.E.S, 'SCOPE', scope):
                root = Path(directory); decoded = saved_pair(root, contigs, calls)
                if change == 'count': decoded['primary_CDS_chain_predictions']['A44_JOINT'] += 1
                elif change == 'GFF_phase':
                    path = root/'A44_JOINT/A44_JOINT.gff3'; lines = path.read_text().splitlines()
                    i = next(i for i, line in enumerate(lines) if '\tCDS\t' in line)
                    fields = lines[i].split('\t'); fields[7] = '1'; lines[i] = '\t'.join(fields); path.write_text('\n'.join(lines)+'\n')
                else:
                    path = root/decoded['scopes']['A44_JOINT'][0]['file']; row = json.loads(path.read_text())
                    row['source_fit_job' if change == 'source' else 'chain_gains'] = '13424143' if change == 'source' else [0.]
                    path.write_text(json.dumps(row)+'\n')
                inputs, reference, comparisons = Mock(return_value=raw), Mock(return_value=context), Mock()
                with self.assertRaises(ValueError): W.evaluate_saved(decoded, root, contigs, inputs, reference, comparisons)
                inputs.assert_not_called(); reference.assert_not_called(); comparisons.assert_not_called()

    def test_both_full_denominators_background_worst_species_and_paired_comparison(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        calls = {('A44_LOCAL', first, '+'): [[(0, 6)]], ('A44_JOINT', first, '+'): [[(0, 6)], [(6, 12)]]}
        with tempfile.TemporaryDirectory() as directory, patch.object(W.V.E.S, 'SCOPE', scope):
            root = Path(directory); decoded = saved_pair(root, contigs, calls)
            predictions, counts = W.reconstruct(decoded, root, contigs)
            current = W.V.measure(predictions['A44_LOCAL'], counts['A44_LOCAL'], context, root, contigs, raw)
            inputs, reference, comparisons = Mock(return_value=raw), Mock(return_value=context), Mock(return_value=comparators(current))
            result = W.evaluate_saved(decoded, root, contigs, inputs, reference, comparisons)
            inputs.assert_called_once(); reference.assert_called_once(); comparisons.assert_called_once()
            methods = result['methods']; table = result['exact_chain_comparison']
            self.assertEqual(set(table), {'A44_LOCAL', 'A44_JOINT', 'A31_A3', 'A31_A23', 'ANNEVO', 'Tiberius', 'Helixer'})
            self.assertEqual(methods['A44_LOCAL']['pooled']['exact_chain'], W.V.E.S.prf(1, 1, 7728))
            self.assertEqual(methods['A44_JOINT']['pooled']['exact_chain'], W.V.E.S.prf(1, 2, 7728))
            self.assertEqual(methods['A44_JOINT']['pooled']['background']['predicted_span_bp'], 6)
            self.assertEqual(table['A44_JOINT']['worst_species_F1'], 0.)
            self.assertEqual(table['A44_JOINT']['macro_F1'], sum(v['exact_chain']['f1'] for v in methods['A44_JOINT']['per_species'].values())/2)
            for arm in W.I.ARMS:
                self.assertEqual(methods[arm]['GFF3'], str(root/arm/(arm+'.gff3')))
                self.assertEqual(methods[arm]['pooled']['exact_chain']['reference'], 7728)
            self.assertIn('paired_JOINT_vs_LOCAL', result); self.assertFalse(result['models_ensembled'])
            old, alternative, native = comparators(current); species = next(iter(methods['A44_LOCAL']['per_species']))
            native['species'][species]['methods']['ANNEVO']['background']['bp'] -= 1
            with self.assertRaises(ValueError): W.compare(methods, old, alternative, native)

    def test_incomplete_explicit_FASTA_precedes_reference_and_full_DNA_mismatch_is_not_filtered(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        with tempfile.TemporaryDirectory() as directory, patch.object(W.V.E.S, 'SCOPE', scope):
            root = Path(directory); decoded = saved_pair(root, contigs, {})
            reference, comparisons = Mock(return_value=context), Mock()
            with self.assertRaisesRegex(ValueError, 'before reference'):
                W.evaluate_saved(decoded, root, contigs, lambda _: {}, reference, comparisons)
            reference.assert_not_called(); comparisons.assert_not_called()
            wrong = dict(raw, **{first: raw[first][:-1]+'A'})
            with self.assertRaisesRegex(ValueError, 'raw FASTA differs'):
                W.evaluate_saved(decoded, root, contigs, lambda _: wrong, reference, comparisons)
            comparisons.assert_not_called(); self.assertEqual(sum(len(v) for v in context[0].values()), 7728)


if __name__ == '__main__': unittest.main()
