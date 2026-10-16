#!/usr/bin/env python3
"""Artificial GFF3/metric composition checks for the new TRAIN diagnostic."""
import json
import os
from pathlib import Path
import tempfile
import time
import evaluate_full_grammar_original_train as E
from src.m28.export import as_transcript, write_gff3

OUT=Path('outputs/M28-FULL-GRAMMAR-TRAIN-DIAGNOSTIC-CHECKS')
C=E.C


def parsed_predictions(directory,seqid,length,txs):
    path=directory/(seqid+'.gff3')
    write_gff3(path,txs,{seqid:length},'ARTIFICIAL_DIAGNOSTIC_CHECK')
    parsed=C.S.E.primary_transcripts(C.S.E.parse_annotation(path,{seqid:length}),complete_only=False)
    expected={t['id']:(t['seqid'],t['strand'],t['CDS']) for t in txs}
    actual={t['id']:(t['seqid'],t['strand'],t['CDS']) for t in parsed}
    assert actual==expected, 'GFF3 one-based inclusive / half-open CDS-phase roundtrip changed'
    return parsed


def main():
    began=time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='m28_train_metric_') as temporary:
        directory=Path(temporary)
        a='artificial_A'; b='artificial_B'
        refs_a=[as_transcript([(10,19)],a,'+','ref_A'),
                as_transcript([(40,49)],a,'+','ref_B')]
        alternate=as_transcript([(40,46),(56,59)],a,'+','ref_B_alt')
        predictions_a=[as_transcript([(10,19)],a,'+','prediction_exact'),
            as_transcript([(40,46),(56,59)],a,'+','prediction_other_isoform'),
            as_transcript([(10,19)],a,'-','prediction_opposite_strand'),
            as_transcript([(90,99)],a,'+','prediction_background')]
        parsed_a=parsed_predictions(directory,a,120,predictions_a)
        # Gene spans deliberately exceed CDS spans: the background denominator
        # must use annotated genes, not reference CDS or model predictions.
        annotation_a={'genes':{'A':{'start':5,'end':30},'B':{'start':35,'end':65}}}
        row_a=C.score(parsed_a,refs_a,refs_a+[alternate],120,annotation_a)
        assert {k:row_a['exact_chain'][k] for k in ('tp','predicted','reference')}==dict(tp=1,predicted=4,reference=2)
        assert row_a['exact_chain']['f1']==1/3
        assert row_a['exact_other_assessable_isoform_not_primary_TP']==1
        assert row_a['duplicate_chain_records']==0
        assert {k:row_a['CDS_base_unstranded_union'][k] for k in ('tp','predicted','reference')}==dict(tp=15,predicted=27,reference=18)
        assert {k:row_a['background'][k] for k in ('bp','predicted_span_bp','predicted_CDS_bp','wholly_background_unique_chains')}==dict(
            bp=65,predicted_span_bp=9,predicted_CDS_bp=9,wholly_background_unique_chains=1)
        refs_b=[as_transcript([(10,19)],b,'+','ref_A')]
        predictions_b=[as_transcript([(10,19)],b,'+','prediction_exact')]
        parsed_b=parsed_predictions(directory,b,80,predictions_b)
        assert C.S.key(parsed_a[0])!=C.S.key(parsed_b[0]), 'Different contigs must retain distinct chain identities'
        row_b=C.score(parsed_b,refs_b,refs_b,80,{'genes':{'A':{'start':5,'end':30}}})
        pooled=E.aggregate([row_a,row_b])
        assert {k:pooled['exact_chain'][k] for k in ('tp','predicted','reference')}==dict(tp=2,predicted=5,reference=3)
        assert pooled['exact_chain']['f1']==0.5
        assert pooled['exact_chain']['f1']!=(row_a['exact_chain']['f1']+row_b['exact_chain']['f1'])/2
        assert {k:pooled['CDS_base_unstranded_union'][k] for k in ('tp','predicted','reference')}==dict(tp=24,predicted=36,reference=27)
        assert pooled['background']['bp']==120 and pooled['background']['span_FPR']==0.075
        assert pooled['exact_other_assessable_isoform_not_primary_TP']==1
    summary=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        scope='artificial_two_contig_GFF3_and_metric_composition_only',
        checks=['GFF3_coordinate_phase_roundtrip_strand_sensitive_exact_chain_other_isoform_and_gene_background',
                'cross_contig_count_pooled_micro_not_mean_F1'],
        artificial_chromosome_A=row_a,artificial_chromosome_B=row_b,artificial_pooled=pooled,
        research_DNA_or_labels_opened=False,checkpoint_loaded=False,actual_native_GLM_forward=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,
        interpretation='evaluation_implementation_check_not_model_accuracy_or_original_TRAIN_result')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(summary,allow_nan=False),flush=True)


if __name__=='__main__': main()
