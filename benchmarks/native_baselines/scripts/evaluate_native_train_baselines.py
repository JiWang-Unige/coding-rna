#!/usr/bin/env python3
"""Original TRAIN native comparisons after all eighteen sequence-only runs finish."""
import json
import os
import time
import run_native_train_baseline as N
import evaluate_full_grammar_original_train as E

C = E.C
ROOT = N.ROOT
OUT = ROOT/'outputs/M28-NATIVE-TRAIN-BASELINE-EVAL'
CONTROL = E.OUT
EXPECTED = {'NC_003070.9':7135,'NC_003071.7':4293,'NC_003076.8':6311,
            'NC_089043.1':1638,'NC_089044.1':1601,'NC_089046.1':1710}


def reference_metadata(row):
    return (row['reference_primary_records'],row['reference_assessable_isoform_records'],
            row['retained_unsupported_reference_strata'],row['background']['bp'],
            row['CDS_base_unstranded_union']['reference'])


def main():
    began = time.perf_counter()
    runs = N.completed_runs()  # No reference or metric file before all native runs finish.
    if (CONTROL/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Original TRAIN comparison control is incomplete')
    control = json.loads((CONTROL/'summary.json').read_text())
    if (not control['completed'] or control['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
            or set(control['chromosomes']) != set(EXPECTED)
            or control['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'):
        raise ValueError('Original TRAIN control reference identity differs')
    report = dict(completed=True,status='native_original_TRAIN_comparison_not_DEV_or_generalization',
        job_id=os.environ['SLURM_JOB_ID'],reference_policy=control['reference_policy'],
        reference_primary_chains=sum(EXPECTED.values()),unsupported_reference_strata_retained=True,
        training_eligibility_does_not_filter_reference=True,reference_used_in_inference=False,
        DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0,input_mode=N.INPUT_MODE,
        control_receipt=str((CONTROL/'summary.json').relative_to(ROOT)),methods={})
    rows_by_method = {method:[] for method in N.METHODS}
    for species,seqid,length in N.SCOPE:
        sequences = C.M.read_allowed(N.input_fasta(seqid),{seqid:length})
        data = ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length})
        with (data/'reference.gff3').open() as handle:
            _,_,raw_cds,_ = C.P.raw_features(handle,seqid)
        refs, isoforms = C.select_reference(coding,raw_cds,sequences)
        if len(refs) != EXPECTED[seqid] or len({C.S.key(t) for t in refs}) != EXPECTED[seqid]:
            raise ValueError('Original TRAIN reference denominator differs; do not filter')
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Reference phase differs; do not exclude references')
        overlaps = C.M.C.overlapping_ids(refs)
        unsupported = C.Counter()
        for tx in refs:
            reasons,_ = C.M.chain_checks(tx,sequences[seqid])
            if tx['id'] in overlaps: reasons.append('same_strand_primary_overlap')
            if len(tx['CDS']) > 128: reasons.append('over_128_exons')
            if tx['CDS'][-1][1]-tx['CDS'][0][0] > C.M.WIDTH: reasons.append('longer_than_native_window')
            unsupported.update(set(reasons))
        for method in N.METHODS:
            folder = N.OUT/method/seqid
            native = runs[(method,seqid)]
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(folder/'prediction.gff3',{seqid:length}),complete_only=False)
            if native['output_check']['feature_counts'].get('CDS',0) > 0 and not txs:
                raise ValueError('Native CDS rows were lost by common parser')
            row = C.score(txs,refs,isoforms,length,annotation)
            row['prediction_nonexact_primary'] = row.pop('reference_nonexact_primary')
            row.update(species=species,seqid=seqid,length=length,reference_primary_records=len(refs),
                reference_assessable_isoform_records=len(isoforms),retained_unsupported_reference_strata=dict(unsupported),
                phase_chain_consistent_prediction_records=sum(C.R7.phase_chain_consistent(t) for t in txs),
                raw_prediction_path=str((folder/'prediction.gff3').relative_to(ROOT)),
                native_run_receipt=str((folder/'summary.json').relative_to(ROOT)),native_runtime_seconds=native['seconds'],
                native_GFF_feature_counts=native['output_check']['feature_counts'],
                baseline_asset_receipt=native['baseline_asset_receipt'])
            if reference_metadata(row) != reference_metadata(control['chromosomes'][seqid]):
                raise ValueError('Reference, retained strata or background differs from original TRAIN control')
            rows_by_method[method].append(row)
        print(json.dumps({'completed_native_TRAIN_scoring':seqid,'reference_primary_chains':len(refs)}),flush=True)
    for method,rows in rows_by_method.items():
        species = sorted({r['species'] for r in rows})
        per_species = {s:E.aggregate([r for r in rows if r['species']==s]) for s in species}
        pooled = E.aggregate(rows)
        f1s = [per_species[s]['exact_chain']['f1'] for s in species]
        macro, worst = sum(f1s)/len(f1s), min(f1s)
        report['methods'][method] = dict(chromosomes={r['seqid']:r for r in rows},per_species=per_species,pooled=pooled,
            macro_species_exact_chain_F1=macro,worst_species_exact_chain_F1=worst,
            native_minus_original4608_F1_pp=dict(
                per_species={s:100*(per_species[s]['exact_chain']['f1']-control['per_species'][s]['exact_chain']['f1']) for s in species},
                pooled=100*(pooled['exact_chain']['f1']-control['pooled']['exact_chain']['f1']),
                macro_species=100*(macro-control['macro_species_exact_chain_F1']),
                worst_species=100*(worst-control['worst_species_exact_chain_F1'])))
    report['evaluation_seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({'completed':True,'reference_primary_chains':report['reference_primary_chains']}),flush=True)


if __name__ == '__main__': main()
