#!/usr/bin/env python3
"""Post-inference diagnostic on original TRAIN only, with the existing CDS ruler."""
import json
import os
from pathlib import Path
import time
import common_ruler as C
import infer_full_grammar_original_train as I

OUT=I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC'


def aggregate(rows):
    result={name:C.S.prf(*(sum(r[name][k] for r in rows) for k in ('tp','predicted','reference')))
            for name in ('exact_chain','CDS_base_unstranded_union')}
    bg={k:sum(r['background'][k] for r in rows) for k in
        ('bp','predicted_span_bp','predicted_CDS_bp','wholly_background_unique_chains')}
    bg['span_FPR']=bg['predicted_span_bp']/bg['bp'] if bg['bp'] else None
    bg['wholly_background_chains_per_Mb']=bg['wholly_background_unique_chains']*1e6/bg['bp'] if bg['bp'] else None
    result.update(background=bg,prediction_records=sum(r['prediction_records'] for r in rows),
        duplicate_chain_records=sum(r['duplicate_chain_records'] for r in rows),
        exact_other_assessable_isoform_not_primary_TP=sum(r['exact_other_assessable_isoform_not_primary_TP'] for r in rows))
    return result


def reconcile_inputs(infer_root=I.OUT, checkpoint=I.CHECKPOINT, identity=I.R.IDENTITY):
    runs={}
    for task,(species,seqid,length) in enumerate(I.SCOPE):
        folder=infer_root/seqid
        if (folder/'STATUS').read_text().strip()!='COMPLETED':
            raise ValueError('All six TRAIN inference tasks must complete before reference access')
        outer=json.loads((folder/'summary.json').read_text())
        raw=json.loads((folder/'prediction/summary.json').read_text())
        if (not outer['inference_complete'] or (outer['task'],outer['species'],outer['seqid'],outer['length'])!=(task,species,seqid,length)
                or outer['checkpoint']!=str(checkpoint) or outer['DEV_used'] or outer['sealed_data_used']
                or outer['reference_used'] or outer['labels_or_feature_cache_used'] or outer['optimizer_steps_this_run']!=0):
            raise ValueError('Original TRAIN-only inference identity differs')
        expected=[[a,min(a+C.M.WIDTH,length)] for a in C.M.starts(length)]
        if (not raw['inference_complete'] or raw['checkpoint']!=str(checkpoint)
                or raw['lengths']!={seqid:length} or raw['reference_used'] or raw['labels_or_feature_cache_used']
                or raw['observation_model_used_at_inference'] or raw['optimizer_steps_this_run']!=0
                or raw['backend']!='exact_exon_budget_skyline'
                or any(raw['configuration'].get(k)!=v for k,v in identity.items())
                or len(raw['scopes'])!=2 or {s['strand'] for s in raw['scopes']}!={'+','-'}):
            raise ValueError('Fixed final-checkpoint whole-contig raw identity differs')
        for scope in raw['scopes']:
            if (scope['seqid']!=seqid or scope['windows']!=expected or scope['reference_used']
                    or not scope['joint_whole_oriented_contig_decode'] or scope['context_policy']!=I.R.CONTEXT
                    or scope['pruning']!='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only'
                    or scope['score_accumulation']!='FP64_per_base_forward_order'):
                raise ValueError('Whole-contig coverage, grammar or scoring differs')
        gff=folder/'prediction/full_grammar.gff3'
        if not gff.is_file() or not gff.stat().st_size: raise ValueError('Raw GFF3 missing')
        runs[seqid]=raw
    return runs


def main(infer_root=I.OUT, checkpoint=I.CHECKPOINT, identity=I.R.IDENTITY, out=OUT):
    began=time.perf_counter()
    runs=reconcile_inputs(infer_root,checkpoint,identity)  # References only after all predictions complete.
    report=dict(completed=True,status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
        job_id=os.environ['SLURM_JOB_ID'],checkpoint=str(checkpoint),
        reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
        reference_used_in_inference=False,training_eligibility_does_not_filter_reference=True,
        unsupported_reference_strata_retained=True,DEV_used=False,sealed_data_used=False,
        optimizer_steps_this_run=0,chromosomes={},per_species={})
    by_species={}
    for species,seqid,length in I.SCOPE:
        folder=infer_root/seqid
        sequences=C.M.read_allowed(folder/'input.fasta',{seqid:length})
        data=I.ROOT/'data/m1_screen'/species
        coding=C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        annotation=C.S.E.parse_annotation(data/'reference.gff3',{seqid:length})
        with (data/'reference.gff3').open() as handle:
            _,_,raw_cds,_=C.P.raw_features(handle,seqid)
        refs,isoforms=C.select_reference(coding,raw_cds,sequences)
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Assessable primary reference phase inconsistent; do not silently filter')
        txs=C.S.E.primary_transcripts(C.S.E.parse_annotation(folder/'prediction/full_grammar.gff3',{seqid:length}),complete_only=False)
        if len(txs)!=runs[seqid]['primary_CDS_chain_predictions']:
            raise ValueError('GFF3 roundtrip prediction count differs')
        for tx in txs:
            failures,_=C.M.chain_checks(tx,sequences[seqid])
            if failures or len(tx['CDS'])>128: raise ValueError('Exported prediction DNA/phase/cap differs')
        for strand in ('+','-'):
            ordered=sorted((t['CDS'][0][0],t['CDS'][-1][1]) for t in txs if t['strand']==strand)
            if any(b>c for (_,b),(c,_) in zip(ordered,ordered[1:])):
                raise ValueError('Exported selected same-strand genes overlap')
        overlaps=C.M.C.overlapping_ids(refs); unsupported=C.Counter()
        for tx in refs:
            failures,_=C.M.chain_checks(tx,sequences[seqid])
            if tx['id'] in overlaps: failures.append('same_strand_primary_overlap')
            if len(tx['CDS'])>128: failures.append('over_128_exons')
            if tx['CDS'][-1][1]-tx['CDS'][0][0]>C.M.WIDTH: failures.append('longer_than_native_window')
            unsupported.update(set(failures))
        row=C.score(txs,refs,isoforms,length,annotation)
        row['prediction_nonexact_primary']=row.pop('reference_nonexact_primary')
        row.update(species=species,seqid=seqid,length=length,reference_primary_records=len(refs),
            reference_assessable_isoform_records=len(isoforms),
            retained_unsupported_reference_strata=dict(unsupported),
            raw_prediction_path=str((folder/'prediction/full_grammar.gff3').relative_to(I.ROOT)),
            phase_checked_prediction_records=len(txs),
            raw_inference_loop_seconds=runs[seqid]['wall_loop_seconds'],
            native_encoding_seconds=sum(s['encoding_seconds'] for s in runs[seqid]['scopes']),
            joint_decoding_seconds=sum(s['decode_seconds'] for s in runs[seqid]['scopes']))
        report['chromosomes'][seqid]=row; by_species.setdefault(species,[]).append(row)
        print(json.dumps(dict(completed_TRAIN_diagnostic=seqid,reference_primary_records=len(refs))),flush=True)
    report['per_species']={species:aggregate(rows) for species,rows in by_species.items()}
    report['pooled']=aggregate(list(report['chromosomes'].values()))
    f1s=[r['exact_chain']['f1'] for r in report['per_species'].values()]
    report['macro_species_exact_chain_F1']=sum(f1s)/len(f1s)
    report['worst_species_exact_chain_F1']=min(f1s)
    report['evaluation_seconds']=time.perf_counter()-began
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True,status=report['status'])),flush=True)


if __name__=='__main__': main()
