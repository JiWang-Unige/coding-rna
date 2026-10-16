"""A43 whole-contig ruler adapter on artificial JSON/GFF/DNA only."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import evaluate_natural_competitor as V
from src.m28.export import as_transcript, write_gff3


def saved_output(root, contigs, calls):
    scopes, txs = [], []
    for ordinal, (seqid, c) in enumerate(contigs.items()):
        length = c['length']
        windows = [[a, min(a+V.I.M.WIDTH, length)] for a in V.I.M.starts(length)]
        for strand, tag in (('+', 'p'), ('-', 'm')):
            chains = calls.get((seqid, strand), [])
            chains = [[list(p) for p in chain] for chain in chains]
            gains = [float(i+1) for i in range(len(chains))]
            converted = [[list(p) for p in V.genomic_chain(c, 0, length, strand)] for c in chains]
            name = 'contig_'+str(ordinal).zfill(4)+'_'+tag+'.json'
            scope = dict(species='user_input', seqid=seqid, strand=strand, ordinal=ordinal,
                windows=windows, file=name, context_policy=V.I.CONTEXT_POLICY,
                primary_CDS_chain_predictions=len(chains), seconds=0.)
            row = dict(variant='A43', score_kind=V.I.SCORE_KIND, species='user_input', seqid=seqid,
                strand=strand, length=length, windows=windows, chains=chains, chain_gains=gains,
                gain=sum(gains), genomic_chains=converted, counts={}, reference_used=False,
                reference_used_in_scoring=False, context_policy=V.I.CONTEXT_POLICY,
                genomic_context_segments=[list(s) for s in V.context_segments(length, windows)])
            (root/name).write_text(json.dumps(row)+'\n'); scopes.append(scope)
            for i, (chain, gain) in enumerate(zip(converted, gains)):
                txs.append(as_transcript(chain, seqid, strand, 'M28_A43_'+str(ordinal)+'_'+tag+'_'+str(i), gain))
    lengths = {q: c['length'] for q, c in contigs.items()}
    write_gff3(root/'A43.gff3', txs, lengths, 'M28_A43')
    decoded = dict(arm='A43', completed=True, inference_complete=True, TRAIN_configuration=V.I.D.CONFIG,
        protocol=V.I.RAW_PROTOCOL, TRAIN_protocol=V.I.D.PROTOCOL,
        source_TRAIN_jobs=dict(fit=V.I.D.FIT_JOB, potentials=V.I.POTENTIAL_JOB, decode=V.I.DECODE_JOB,
            readout=V.I.READOUT_JOB, CPU=V.I.CHECK_JOB), checkpoint=V.I.D.CONFIG['source_checkpoint'],
        optimizer_steps_source=4608, model_parameters=366191, feature_revision=V.I.NATIVE.F.REVISION,
        feature_spec=V.I.D.R.C.R.FEATURE_SPEC, input_mode='explicit_raw_FASTA',
        input_fasta='/synthetic/not_opened.fa', lengths=lengths, scopes=scopes,
        context_policy=V.I.CONTEXT_POLICY, score_kind=V.I.SCORE_KIND, max_exons=128, threshold=0., null=0.,
        optimizer_steps_this_run=0, new_parameters=0, primary_CDS_chain_predictions=len(txs),
        **{k: False for k in ('reference_used', 'labels_or_feature_cache_used', 'generator_model_loaded',
            'query_or_beam_decoding', 'local_selected_only_owner_assembly', 'pretrained_encoder_updated')})
    return decoded, txs


def full_context():
    """7728 artificial reference chains; NOT actual DEV records or measurements."""
    contigs, refs, isoforms, parents, annotations, sequences, lengths, scope = {}, {}, {}, {}, {}, {}, {}, []
    for index, (species, number) in enumerate(V.A.REF_COUNTS.items()):
        seqid, length = 'toy_'+str(index), number*20+100
        contigs[seqid] = dict(species=species, length=length); scope.append((species, seqid, length, number))
        rows = [as_transcript([(i*20, i*20+6)], seqid, '+', 'truth_'+str(index)+'_'+str(i)) for i in range(number)]
        # One reference exceeds the decoder window width; keep it in the full denominator.
        rows[1] = as_transcript([(20, 20+V.I.M.WIDTH+6)], seqid, '+', 'truth_long_'+str(index))
        refs[species] = rows; isoforms[species] = rows; parents[species] = rows
        annotations[species] = dict(genes={t['gene_id']: dict(start=t['CDS'][0][0], end=t['CDS'][-1][1]) for t in rows})
        sequences[species] = {seqid: 'ATGTAAATGTAA'+'C'*(length-12)}; lengths[seqid] = length
    context = refs, isoforms, parents, annotations, sequences, lengths
    raw = {q: sequences[c['species']][q] for q, c in contigs.items()}
    return contigs, context, raw, scope


def comparators(current):
    common = dict(reference_primary_chains=7728, test_setaria_used=False)
    previous = dict(common, methods=dict(A31_A3=deepcopy(current)))
    alternative = dict(common, methods=dict(A31_A23=deepcopy(current)))
    native = dict(common, completed=True, unsupported_reference_strata_retained=True,
        species={s: dict(methods={n: deepcopy(v) for n in ('ANNEVO', 'Tiberius', 'Helixer')})
            for s, v in current['per_species'].items()},
        pooled={n: deepcopy(current['pooled']) for n in ('ANNEVO', 'Tiberius', 'Helixer')},
        source_limitations=dict(synthetic=True, not_actual_native_results=True))
    return previous, alternative, native


class NaturalCompetitorEvaluationTests(unittest.TestCase):
    def test_actual_single_model_receipt_not_old_pool_or_implicit_species(self):
        contigs = dict(toy=dict(species='synthetic', length=40))
        with tempfile.TemporaryDirectory() as directory:
            decoded, _ = saved_output(Path(directory), contigs, {})
            V.validate_complete(decoded, contigs)
            for key, value in (('checkpoint', 'early.pt'), ('model_parameters', 649457),
                    ('input_mode', 'fixed_DEV_raw_FASTA'), ('generator_model_loaded', True),
                    ('local_selected_only_owner_assembly', True), ('source_TRAIN_jobs', {'fit': 'other'}),
                    ('protocol', 'old_protocol'), ('threshold', .1), ('reference_used', True)):
                with self.subTest(key=key), self.assertRaises(ValueError):
                    V.validate_complete(dict(decoded, **{key: value}), contigs)
            with self.assertRaisesRegex(ValueError, 'mapping'):
                V.validate_complete(decoded, {})

    def test_complete_both_strand_original_grid_ordinal_and_explicit_FASTA_species(self):
        contigs = dict(toy=dict(species='synthetic', length=36864))
        with tempfile.TemporaryDirectory() as directory:
            decoded, _ = saved_output(Path(directory), contigs, {})
            for change in ('missing', 'duplicate', 'grid', 'ordinal', 'species', 'context'):
                bad = deepcopy(decoded)
                if change == 'missing': bad['scopes'].pop()
                elif change == 'duplicate': bad['scopes'][1] = deepcopy(bad['scopes'][0])
                elif change == 'grid': bad['scopes'][0]['windows'][1][0] += 1
                elif change == 'ordinal': bad['scopes'][1]['ordinal'] = 1
                elif change == 'species': bad['scopes'][0]['species'] = 'arabidopsis_thaliana'
                else: bad['scopes'][0]['context_policy'] = 'old_chain_owner'
                with self.subTest(change=change), self.assertRaises(ValueError): V.validate_complete(bad, contigs)

    def test_minus_chain_order_gain_pairing_and_GFF_phase_before_reference(self):
        contigs = dict(toy=dict(species='synthetic', length=40))
        calls = {('toy', '-'): [[(0, 4), (9, 14)], [(20, 26)]]}
        with tempfile.TemporaryDirectory() as directory, patch.object(V.A, 'reference_context') as reference,\
                patch('src.m28.assembly.assemble_predictions', side_effect=AssertionError('No old owner assembly')):
            root = Path(directory); decoded, expected = saved_output(root, contigs, calls)
            predictions, counts = V.reconstruct(decoded, root, contigs)
            self.assertEqual(predictions['synthetic'], expected)
            self.assertEqual(predictions['synthetic'][0]['CDS'], [(26, 31, '2'), (36, 40, '0')])
            self.assertEqual([t['score'] for t in predictions['synthetic']], [1., 2.])
            self.assertGreater(predictions['synthetic'][0]['CDS'][0][0], predictions['synthetic'][1]['CDS'][0][0])
            self.assertEqual(counts['toy|-']['selected_chains'], 2); reference.assert_not_called()

    def test_joint_geometry_context_nonoverlap_cap_positive_gains_and_totals(self):
        contigs = dict(toy=dict(species='synthetic', length=40))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); decoded, _ = saved_output(root, contigs, {('toy', '+'): [[(0, 6)], [(20, 26)]]})
            scope = decoded['scopes'][0]; row = json.loads((root/scope['file']).read_text())
            V.validate_row(row, scope, 40)
            for key, value in (('chain_gains', [0., 2.]), ('chain_gains', [float('nan'), 2.]), ('gain', 9.),
                    ('chains', [[[0, 6]], [[5, 11]]]), ('chains', [[[0, 4], [8, 13]], [[20, 26]]]),
                    ('chains', [[[0., 6]], [[20, 26]]]), ('chains', [[[0, 6]]*129, [[20, 26]]]),
                    ('genomic_chains', list(reversed(row['genomic_chains']))),
                    ('genomic_context_segments', [[0, 0, 39]]), ('reference_used_in_scoring', True)):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    V.validate_row(dict(row, **{key: value}), scope, 40)

    def test_GFF_coordinates_phase_and_summary_count_changes_are_not_accepted(self):
        contigs = dict(toy=dict(species='synthetic', length=40))
        for change in ('coordinate', 'phase', 'count'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); decoded, _ = saved_output(root, contigs, {('toy', '-'): [[(0, 4), (9, 14)]]})
                if change == 'count': decoded['primary_CDS_chain_predictions'] += 1
                else:
                    lines = (root/'A43.gff3').read_text().splitlines()
                    index = next(i for i, line in enumerate(lines) if '\tCDS\t' in line)
                    fields = lines[index].split('\t'); fields[3 if change == 'coordinate' else 7] = '28' if change == 'coordinate' else '1'
                    lines[index] = '\t'.join(fields); (root/'A43.gff3').write_text('\n'.join(lines)+'\n')
                with self.assertRaises(ValueError): V.reconstruct(decoded, root, contigs)

    def test_cross_window_long_chain_and_zero_predictions_preserve_complete_scopes(self):
        contigs = dict(long=dict(species='synthetic', length=36864), empty=dict(species='other', length=40))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); decoded, _ = saved_output(root, contigs, {('long', '+'): [[(3000, 33000)]]})
            predictions, counts = V.reconstruct(decoded, root, contigs)
            self.assertEqual(predictions['synthetic'][0]['CDS'], [(3000, 33000, '0')])
            self.assertEqual(predictions['other'], [])
            self.assertEqual(len(counts), 4); self.assertEqual(counts['empty|+']['selected_gain'], 0.)

    def test_full_7728_synthetic_denominators_background_and_native_worst_species(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        with tempfile.TemporaryDirectory() as directory, patch.object(V.E.S, 'SCOPE', scope):
            root = Path(directory); decoded, _ = saved_output(root, contigs, {(first, '+'): [[(0, 6)], [(6, 12)]]})
            predictions, counts = V.reconstruct(decoded, root, contigs)
            current = V.measure(predictions, counts, context, root, contigs, raw)
            self.assertEqual(current['pooled']['exact_chain'], V.E.S.prf(1, 2, 7728))
            self.assertEqual(current['pooled']['background']['predicted_span_bp'], 6)
            self.assertEqual(sum(len(v) for v in context[0].values()), 7728)
            self.assertTrue(any(t['CDS'][-1][1]-t['CDS'][0][0] > V.I.M.WIDTH for ts in context[0].values() for t in ts))
            self.assertNotIn('actual_free_primary_coverage', current['pooled'])
            compared = V.compare(current, *comparators(current)); table = compared['exact_chain_comparison']
            self.assertEqual(set(table), {'A43', 'A31_A3', 'A31_A23', 'ANNEVO', 'Tiberius', 'Helixer'})
            self.assertEqual(table['A43']['worst_species_F1'], 0.)
            self.assertEqual(table['A43']['macro_F1'], sum(v['exact_chain']['f1'] for v in current['per_species'].values())/2)
            for change in ('reference', 'background', 'unsupported'):
                old, alternative, native = comparators(current); species = next(iter(current['per_species']))
                if change == 'unsupported': native['unsupported_reference_strata_retained'] = False
                else: native['species'][species]['methods']['ANNEVO']['exact_chain' if change == 'reference' else 'background']['reference' if change == 'reference' else 'bp'] -= 1
                with self.subTest(change=change), self.assertRaises(ValueError): V.compare(current, old, alternative, native)

    def test_DNA_mismatch_or_illegal_prediction_stops_without_truth_filtering(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        with tempfile.TemporaryDirectory() as directory, patch.object(V.E.S, 'SCOPE', scope):
            root = Path(directory); decoded, _ = saved_output(root, contigs, {(first, '+'): [[(0, 6)]]})
            predictions, counts = V.reconstruct(decoded, root, contigs)
            wrong = dict(raw, **{first: raw[first][:-1]+'A'})
            with self.assertRaisesRegex(ValueError, 'raw FASTA differs'):
                V.measure(predictions, counts, context, root, contigs, wrong)
            bad = deepcopy(predictions); species = contigs[first]['species']; bad[species][0]['CDS'] = [(1, 7, '0')]
            with self.assertRaisesRegex(ValueError, 'grammar'):
                V.measure(bad, counts, context, root, contigs, raw)
            self.assertEqual(sum(len(v) for v in context[0].values()), 7728)

    def test_saved_closure_and_complete_input_precede_reference_and_comparison_access(self):
        contigs, context, raw, scope = full_context(); first = next(iter(contigs))
        with tempfile.TemporaryDirectory() as directory, patch.object(V.E.S, 'SCOPE', scope):
            root = Path(directory); decoded, _ = saved_output(root, contigs, {(first, '+'): [[(0, 6)]]})
            input_loader, reference_loader, comparison_loader = Mock(return_value=raw), Mock(return_value=context), Mock()
            with self.assertRaises(ValueError):
                V.evaluate_saved(dict(decoded, primary_CDS_chain_predictions=2), root, contigs,
                    input_loader, reference_loader, comparison_loader)
            input_loader.assert_not_called(); reference_loader.assert_not_called(); comparison_loader.assert_not_called()
            with self.assertRaises(ValueError):
                V.evaluate_saved(decoded, root, {first: contigs[first]}, input_loader, reference_loader, comparison_loader)
            reference_loader.assert_not_called()
            with self.assertRaisesRegex(ValueError, 'before reference'):
                V.evaluate_saved(decoded, root, contigs, Mock(return_value={}), reference_loader, comparison_loader)
            reference_loader.assert_not_called(); comparison_loader.assert_not_called()
            pred, counts = V.reconstruct(decoded, root, contigs)
            current = V.measure(pred, counts, context, root, contigs, raw)
            comparison_loader.return_value = comparators(current)
            report = V.evaluate_saved(decoded, root, contigs, input_loader, reference_loader, comparison_loader)
            input_loader.assert_called_once_with('/synthetic/not_opened.fa'); reference_loader.assert_called_once()
            comparison_loader.assert_called_once()
            self.assertFalse(report['finite_pool_owner_coverage_applicable'])
            self.assertTrue(report['reported_raw_FASTA_matches_complete_reference_contig_DNA'])
            self.assertTrue(report['raw_predictions_and_GFF_reconciled_before_reference_access'])
            for key in ('independent_neural_forward_replay', 'independent_saved_emission_score_replay',
                    'independent_full_grammar_DP_replay'): self.assertFalse(report[key])


if __name__ == '__main__':
    unittest.main()
