#!/usr/bin/env python3
"""One common CDS-assessable DEV ruler; cached predictions are never altered."""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import reference_pair as R
import semantic_reference_pair as R7
import manifest as M
P, S = R.P, R.S


def select_reference(coding, raw_cds, sequences):
    eligible = defaultdict(list)
    for tx in coding['transcripts'].values():
        if tx['CDS'] and not any(S.E.is_partial(a) for a in raw_cds.get(tx['id'], [])):
            if P.compatible(tx, S.V.transcript_sequence(tx, sequences).upper()):
                eligible[tx['gene_id']].append(tx)
    return R.pick_primary(eligible), [t for ts in eligible.values() for t in ts]


def score(txs, refs, all_isoforms, length, annotation):
    # Unique complete CDS chains, not transcript IDs; exact phase is diagnosed separately.
    pk, rk, ik = ({S.key(t) for t in rows} for rows in (txs, refs, all_isoforms))
    ps, pc = P.coverage_masks(txs, length)
    _, rc = P.coverage_masks(refs, length)
    genes = np.zeros(length, bool)
    for g in annotation['genes'].values():
        genes[g['start']:g['end']] = True
    bg = ~genes
    wholly_background = sum(not genes[t['CDS'][0][0]:t['CDS'][-1][1]].any()
                            for t in {S.key(t):t for t in txs}.values())
    bg_bp = int(bg.sum())
    return {
        'exact_chain': S.prf(len(pk & rk), len(pk), len(rk)),
        'CDS_base_unstranded_union': S.prf(int((pc & rc).sum()), int(pc.sum()), int(rc.sum())),
        'exact_other_assessable_isoform_not_primary_TP': len((pk & ik)-rk),
        'reference_nonexact_primary': len(pk-rk),
        'prediction_records': len(txs), 'duplicate_chain_records': len(txs)-len(pk),
        'background': {
            'definition': 'outside_all_annotated_gene_feature_spans_either_strand',
            'bp': bg_bp, 'predicted_span_bp': int((ps & bg).sum()),
            'predicted_CDS_bp': int((pc & bg).sum()),
            'span_FPR': float((ps & bg).sum()/bg_bp) if bg_bp else None,
            'wholly_background_unique_chains': wholly_background,
            'wholly_background_chains_per_Mb': wholly_background*1e6/bg_bp if bg_bp else None}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    if (out/'summary.json').exists():
        raise FileExistsError('Completed result exists')
    old = json.loads((ROOT/'outputs/M26-SAME-SCOPE-MECHANISM-R1/summary.json').read_text())
    result = {'status':'development_rescore_not_independent_test',
              'policy':'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
              'input_condition':'saved_sequence_only_predictions_no_new_tool_runs',
              'test_access':False, 'setaria_access':False,
              'training_eligibility_does_not_filter_reference':True, 'species':{}, 'pooled':{}}
    pooled = {m:[] for m in S.METHODS}
    for species, seqid, length, expected_old in S.SCOPE:
        path = ROOT/'data/m1_screen'/species
        sequences = M.read_allowed(path/'genome.fa', {seqid:length})
        coding = S.E.parse_annotation(path/'reference.gff3', {seqid:length}, protein_coding_only=True)
        annotation = S.E.parse_annotation(path/'reference.gff3', {seqid:length})
        with (path/'reference.gff3').open() as f:
            _, _, raw_cds, _ = P.raw_features(f, seqid)
        refs, isoforms = select_reference(coding, raw_cds, sequences)
        if any(not R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('CDS-assessable primary phase mismatch: do not silently drop references')
        parent = S.E.primary_transcripts(coding)
        assert len(parent) == expected_old
        if species == 'arabidopsis_thaliana':
            assert len(refs) == 5437, 'M26 established reference policy did not replay'
        r2 = {t['gene_id']:t for t in S.E.primary_transcripts(coding, complete_only=False)}
        changes = [t['gene_id'] for t in refs if S.key(t) != S.key(r2[t['gene_id']])]
        overlaps = M.C.overlapping_ids(refs)
        failures = Counter()
        with (out/(species+'.references.jsonl')).open('x') as f:
            for t in refs:
                unsupported, motifs = M.chain_checks(t, sequences[seqid])
                if t['id'] in overlaps: unsupported.append('same_strand_primary_overlap')
                if t['CDS'][-1][1]-t['CDS'][0][0] > M.WIDTH:
                    unsupported.append('longer_than_window')
                failures.update(set(unsupported))
                f.write(json.dumps({'species':species, 'id':t['id'], 'gene_id':t['gene_id'],
                                    'seqid':seqid, 'strand':t['strand'], 'CDS':t['CDS'],
                                    'M28_unsupported_strata_not_excluded':sorted(set(unsupported)),
                                    'splice_motifs':motifs})+'\n')
        data = {'seqid':seqid, 'length':length, 'coding_loci':len(coding['genes']),
                'CDS_assessable_primary_records':len(refs),
                'CDS_assessable_unique_primary_chains':len({S.key(t) for t in refs}),
                'CDS_assessable_isoform_records':len(isoforms),
                'historical_parent_primary_records':len(parent),
                'primary_chain_changes_vs_R2_longest_all':len(changes),
                'primary_change_gene_ids':changes,
                'retained_unsupported_reference_strata':dict(failures), 'methods':{}}
        for method, template in S.METHODS.items():
            predpath = ROOT/template.format(species=species)
            txs = S.E.primary_transcripts(S.E.parse_annotation(predpath, {seqid:length}), complete_only=False)
            pk, hk = {S.key(t) for t in txs}, {S.key(t) for t in parent}
            replay = S.prf(len(pk & hk), len(pk), len(hk))
            assert replay == old['methods'][method]['per_species'][species]['exact_chain'], (method,species,'old replay')
            measured = score(txs, refs, isoforms, length, annotation)
            measured['historical_parent_exact_chain_replayed'] = replay
            measured['prediction_path'] = str(predpath.relative_to(ROOT))
            data['methods'][method] = measured
            pooled[method].append(measured)
        result['species'][species] = data
    for method, rows in pooled.items():
        q = {'exact_chain':S.prf(*(sum(r['exact_chain'][k] for r in rows) for k in ('tp','predicted','reference'))),
             'CDS_base_unstranded_union':S.prf(*(sum(r['CDS_base_unstranded_union'][k] for r in rows) for k in ('tp','predicted','reference')))}
        bg = {k:sum(r['background'][k] for r in rows) for k in ('bp','predicted_span_bp','predicted_CDS_bp','wholly_background_unique_chains')}
        bg['span_FPR'] = bg['predicted_span_bp']/bg['bp']
        bg['wholly_background_chains_per_Mb'] = bg['wholly_background_unique_chains']*1e6/bg['bp']
        q['background'] = bg
        q['exact_other_assessable_isoform_not_primary_TP'] = sum(r['exact_other_assessable_isoform_not_primary_TP'] for r in rows)
        result['pooled'][method] = q
    with (out/'summary.json').open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
