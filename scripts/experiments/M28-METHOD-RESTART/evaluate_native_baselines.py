#!/usr/bin/env python3
"""Score the six completed new native predictions with the frozen 7728-chain DEV ruler."""
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common_ruler as C
from run_native_baseline import METHODS,SCOPE,ASSET_PREP,OUT as PRED
OUT=ROOT/'outputs/M28-NATIVE-BASELINE-EVAL'
EXPECTED={'arabidopsis_thaliana':5437,'oryza_sativa':2291}

def main():
    began=time.perf_counter()
    runs={}
    # All inference contracts are reconciled before opening a reference.
    for method in METHODS:
        for species,seqid,length in SCOPE:
            folder=PRED/method/species
            row=json.loads((folder/'summary.json').read_text())
            if (folder/'STATUS').read_text().strip()!='COMPLETED':
                raise ValueError('Incomplete native run '+str(folder))
            if not row['completed'] or (row['method'],row['species'],row['seqid'],row['length'])!=(method,species,seqid,length):
                raise ValueError('Native run identity differs')
            if row['input_mode']!='fixed_DEV_raw_FASTA' or row['optimizer_steps']!=0 or any(
                row[k] for k in ('reference_used','labels_or_feature_cache_used','test_setaria_used')):
                raise ValueError('Native sequence-only boundary differs')
            if row['raw_GFF']!='prediction.gff3' or not (folder/'prediction.gff3').is_file():
                raise ValueError('Native raw GFF missing')
            runs[(method,species)]=row
    report={'completed':True,'status':'native_comparators_two_DEV_not_independent_generalization',
            'reference_primary_chains':sum(EXPECTED.values()),'reference_policy':'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
            'input_condition':'six_fixed_version_native_raw_FASTA_runs_no_external_evidence',
            'unsupported_reference_strata_retained':True,'test_setaria_used':False,
            'species':{},'pooled':{},
            'source_limitations':'pretrained species exposure is method-dependent or unknown; Tiberius author metadata includes rice; no clean universal zero-shot claim',
            'historical_comparison_limit':'old cached predictions retain separate records; version, weight, input set and execution environment differ, so differences are not isolated version effects'}
    pooled={m:[] for m in METHODS}
    for species,seqid,length in SCOPE:
        fasta=ASSET_PREP/'inputs'/(species+'.fa')
        sequence=C.M.read_allowed(fasta,{seqid:length})
        data=ROOT/'data/m1_screen'/species
        coding=C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        annotation=C.S.E.parse_annotation(data/'reference.gff3',{seqid:length})
        with (data/'reference.gff3').open() as handle:
            _,_,raw_cds,_=C.P.raw_features(handle,seqid)
        refs,isoforms=C.select_reference(coding,raw_cds,sequence)
        if len(refs)!=EXPECTED[species] or len({C.S.key(t) for t in refs})!=EXPECTED[species]:
            raise ValueError('Frozen reference denominator differs')
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Frozen reference phase differs; do not exclude affected references')
        overlaps=C.M.C.overlapping_ids(refs)
        unsupported=C.Counter()
        for tx in refs:
            reasons,_=C.M.chain_checks(tx,sequence[seqid])
            if tx['id'] in overlaps: reasons.append('same_strand_primary_overlap')
            if tx['CDS'][-1][1]-tx['CDS'][0][0]>C.M.WIDTH: reasons.append('longer_than_window')
            unsupported.update(set(reasons))
        item={'seqid':seqid,'length':length,'reference_primary_chains':len(refs),
              'retained_unsupported_reference_strata':dict(unsupported),'methods':{}}
        for method in METHODS:
            folder=PRED/method/species
            native=runs[(method,species)]
            txs=C.S.E.primary_transcripts(C.S.E.parse_annotation(folder/'prediction.gff3',{seqid:length}),complete_only=False)
            if native['output_check']['feature_counts'].get('CDS',0)>0 and not txs:
                raise ValueError('Native CDS rows not represented by common parser: '+method+'/'+species)
            measured=C.score(txs,refs,isoforms,length,annotation)
            measured['prediction_nonexact_primary']=measured.pop('reference_nonexact_primary')
            measured['phase_chain_consistent_prediction_records']=sum(C.R7.phase_chain_consistent(t) for t in txs)
            measured['raw_prediction_path']=str((folder/'prediction.gff3').relative_to(ROOT))
            measured['native_run_receipt']=str((folder/'summary.json').relative_to(ROOT))
            measured['native_GFF_feature_counts']=native['output_check']['feature_counts']
            measured['native_runtime_seconds']=native['seconds']
            measured['baseline_asset_receipt']=native['baseline_asset_receipt']
            item['methods'][method]=measured
            pooled[method].append(measured)
        report['species'][species]=item
    for method,rows in pooled.items():
        result={name:C.S.prf(*(sum(r[name][k] for r in rows) for k in ('tp','predicted','reference')))
                for name in ('exact_chain','CDS_base_unstranded_union')}
        result['macro_exact_chain_F1']=sum(r['exact_chain']['f1'] for r in rows)/len(rows)
        bg={k:sum(r['background'][k] for r in rows) for k in
            ('bp','predicted_span_bp','predicted_CDS_bp','wholly_background_unique_chains')}
        bg['span_FPR']=bg['predicted_span_bp']/bg['bp']
        bg['wholly_background_chains_per_Mb']=bg['wholly_background_unique_chains']*1e6/bg['bp']
        result['background']=bg
        report['pooled'][method]=result
    report['evaluation_seconds']=time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({'completed':True,'reference_primary_chains':report['reference_primary_chains'],'methods':list(METHODS),
                      'seconds':report['evaluation_seconds']}),flush=True)

if __name__=='__main__': main()
