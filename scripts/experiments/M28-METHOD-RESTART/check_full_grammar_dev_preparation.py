#!/usr/bin/env python3
"""Artificial-only checks for prepared DEV adapters; never execute research DEV."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
from unittest.mock import patch
import evaluate_full_grammar_dev_baselines as E

OUT = E.I.ROOT/'outputs/M28-FULL-GRAMMAR-DEV-PREPARATION-CHECKS'


def rejects(function,*args):
    try:
        function(*args)
    except ValueError:
        return
    raise AssertionError('Expected the concrete mismatch to be rejected')


def main():
    began = time.perf_counter(); checks = []
    I,R,C = E.I,E.I.R,E.C
    assert I.SCOPE == E.N.SCOPE
    assert E.EXPECTED == {'arabidopsis_thaliana':5437,'oryza_sativa':2291}
    assert I.CHECKPOINT.name == 'step_004608.pt'
    species,seqid,length = 'artificial_species','artificial_dev_contig',40
    outer = dict(inference_complete=True,task=0,species=species,seqid=seqid,length=length,
        checkpoint=str(I.CHECKPOINT),raw_scope=I.RAW_SCOPE,DEV_used=True,sealed_data_used=False,
        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0)
    scope = dict(seqid=seqid,windows=[[a,min(a+C.M.WIDTH,length)] for a in C.M.starts(length)],
        reference_used=False,joint_whole_oriented_contig_decode=True,context_policy=R.CONTEXT,
        pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
        score_accumulation='FP64_per_base_forward_order',
        tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score')
    raw = dict(inference_complete=True,checkpoint=str(I.CHECKPOINT),protocol=R.PROTOCOL,
        lengths={seqid:length},configuration=dict(R.IDENTITY),reference_used=False,
        labels_or_feature_cache_used=False,observation_model_used_at_inference=False,
        optimizer_steps_this_run=0,backend='exact_exon_budget_skyline',
        scopes=[dict(scope,strand='+'),dict(scope,strand='-')])
    E.check_candidate(outer,raw,0,species,seqid,length)
    for field,value in (('checkpoint','step_003072.pt'),('reference_used',True),('sealed_data_used',True)):
        changed = dict(outer); changed[field] = value
        rejects(E.check_candidate,changed,raw,0,species,seqid,length)
    for mutation in ('epsilon','strand','coverage'):
        changed = copy.deepcopy(raw)
        if mutation == 'epsilon': changed['configuration']['epsilon'] = .2
        if mutation == 'strand': changed['scopes'][1]['strand'] = '+'
        if mutation == 'coverage': changed['scopes'][0]['windows'][0][1] -= 1
        rejects(E.check_candidate,outer,changed,0,species,seqid,length)
    checks.append('candidate_final_checkpoint_label_free_fixed_scores_both_strands_and_complete_grid')
    assets = {
        'ANNEVO':dict(source_commit='42c920f70f52f96c7c06c665a639ec47d9304fff',version='2.3.3',
            lineage='Magnoliopsida',weight='assets/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt'),
        'Tiberius':dict(native_metadata=dict(source_commit='871d4979568ae09e56de7eca5d7489046ee833ea',version='2.0.8'),
            base_container_actual_metadata=dict(version='2.0.7'),configuration_commit='871d4979568ae09e56de7eca5d7489046ee833ea',
            entry_point='python3 -m tiberius.main',weights_directory=str(
                (E.N.ASSET_PREP/'assets/Tiberius-model-weights/angiosperms').relative_to(I.ROOT))),
        'Helixer':dict(native_metadata=dict(source_commit='d17bb496b8842542590644b2437798711c0125d5',version='0.3.7'),
            model_source_commit='d17bb496b8842542590644b2437798711c0125d5',
            selected_model=dict(model_file_name='land_plant_v0.3_a_0080.h5',priority='80'))}
    for method,asset in assets.items():
        row = dict(completed=True,method=method,species=species,seqid=seqid,length=length,
            input_mode='fixed_DEV_raw_FASTA',raw_GFF='prediction.gff3',
            case_policy='uppercase_same_as_M28_raw_caller',native_inference_executed=True,
            optimizer_steps=0,reference_used=False,labels_or_feature_cache_used=False,test_setaria_used=False,
            GPU_runtime=dict(visible_gpus=1),baseline_asset_receipt=asset)
        E.check_native_receipt(row,method,species,seqid,length)
        changed = copy.deepcopy(row)
        if method == 'ANNEVO': changed['baseline_asset_receipt']['lineage'] = 'Monocots'
        if method == 'Tiberius': changed['baseline_asset_receipt']['native_metadata']['version'] = '2.0.7'
        if method == 'Helixer': changed['baseline_asset_receipt']['selected_model']['model_file_name'] = 'land_plant_v0.3_m_0100.h5'
        rejects(E.check_native_receipt,changed,method,species,seqid,length)
        rejects(E.check_native_receipt,row,method,'other_species',seqid,length)
    checks.append('pinned_native_receipt_species_versions_lineage_and_weights_not_alternative_assets')
    rows = [[key,'COMPLETED','0:0'] for key in ('17','17.batch','17.extern')]
    E.check_terminal_rows(rows,{'17'})
    rejects(E.check_terminal_rows,rows[:2],{'17'})
    rejects(E.check_terminal_rows,[rows[0],['17.batch','TIMEOUT','0:0'],rows[2]],{'17'})
    rejects(E.check_terminal_rows,[rows[0],['17.batch','COMPLETED','0:1'],rows[2]],{'17'})
    with patch.object(E,'reconcile_inputs',side_effect=ValueError('synthetic incomplete input')), \
         patch.object(C,'select_reference') as select:
        rejects(E.main)
        select.assert_not_called()
    checks.append('complete_successful_Slurm_steps_and_all_inputs_before_reference_selection')
    from src.m28.export import as_transcript
    single = as_transcript([(0,9)],seqid,'+','one')
    multiple = as_transcript([(12,18),(30,36)],seqid,'+','two')
    opposite = as_transcript([(12,18),(30,36)],seqid,'-','opposite')
    measured = E.structure_metrics([single,opposite],[single,multiple])
    assert measured['single_CDS_segment'] == C.S.prf(1,1,1)
    assert measured['multiple_CDS_segments'] == C.S.prf(0,1,1)
    checks.append('structure_strata_keep_strand_and_all_reference_keys')
    with tempfile.TemporaryDirectory(prefix='m28_dev_adapter_artificial_') as directory:
        folder = Path(directory)/seqid; folder.mkdir()
        artificial = 'acgtRynnacgt'
        def read_only_fake(path,allowed):
            assert allowed == {seqid:len(artificial)}
            assert path == I.ROOT/'data/m1_screen'/species/'genome.fa'
            return {seqid:artificial}
        def fake_caller(fasta,checkpoint,out,limit):
            assert fasta.read_text() == '>'+seqid+'\n'+artificial.upper()+'\n'
            assert checkpoint == I.CHECKPOINT and out == folder/'prediction' and limit == .125
        with patch.object(I,'OUT',Path(directory)),patch.object(I,'SCOPE',((species,seqid,len(artificial)),)), \
             patch.object(R.M,'read_allowed',side_effect=read_only_fake) as source, \
             patch.object(R,'main',side_effect=fake_caller) as caller:
            I.main(0)
            source.assert_called_once(); caller.assert_called_once()
            receipt = json.loads((folder/'summary.json').read_text())
            assert receipt['DEV_used'] and not receipt['sealed_data_used'] and not receipt['reference_used']
            rejects(I.main,1)
    checks.append('artificial_driver_preserves_original_uppercase_IUPAC_and_uses_unchanged_raw_caller')
    result = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],checks=checks,
        scope='artificial_DEV_adapter_preparation_only',DEV_used=False,sealed_data_used=False,
        research_DNA_or_labels_opened=False,historical_DEV_receipts_or_predictions_opened=False,
        checkpoint_loaded=False,actual_native_GLM_forward=False,optimizer_steps_this_run=0,
        interpretation='synthetic_preparation_wiring_not_actual_DEV_readiness_or_model_performance',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
