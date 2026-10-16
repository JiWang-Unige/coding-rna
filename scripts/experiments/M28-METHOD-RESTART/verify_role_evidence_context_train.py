#!/usr/bin/env python3
"""A40 actual saved student/source exposure, role inputs, components and gates."""
import argparse
from collections import Counter
import math
import time
import torch
import readout_role_evidence_context_train as D
from verify_context_chain_train import verify_components, expected_optimizer_steps
from src.m28.local_path_score import local_path_scores


@torch.inference_mode()
def verify_fit(fit_job):
    fitted=D.R.completed(D.F.OUT,fit_job); D.validate_fit(fitted)
    draws,_,_,records,_,_,_=D.F.A.A.population()
    training=D.R.C.R.rows(D.F.OUT/'training.jsonl')
    if len(training)!=4608 or len(fitted['epochs'])!=3:
        raise ValueError('Complete actual original three-pass fit required')
    expected_counts=Counter()
    for draw in draws: expected_counts.update(records[draw['window_id']]['counts'])
    for i,row in enumerate(training):
        epoch,ordinal=divmod(i,1536); draw=draws[ordinal]; rec=records[draw['window_id']]
        fields=('weighted_loss','unweighted_mean_known_BCE','unweighted_mean_known_set_NLL',
                'unweighted_fixed_composite_loss','log_partition','gold_gain','gradient_norm_before_clip')
        if (row['step']!=i+1 or row['epoch']!=epoch+1 or row['draw_ordinal']!=ordinal or
                row['window_id']!=draw['window_id'] or row['species']!=draw['species'] or
                row['importance_weight']!=draw['importance_weight_uniform_within_stratum'] or
                any(row[k]!=v for k,v in rec['counts'].items()) or
                any(not math.isfinite(row[k]) for k in fields) or
                row['frozen_local_source_gradients_absent'] is not True or
                (rec['counts']['zero_known_windows'] and any(row[k]!=0 for k in fields))):
            raise ValueError('Actual original draw/weight/counts/finite protected update record differs')
        combined=(row['unweighted_mean_known_BCE']+row['unweighted_mean_known_set_NLL'])/2
        if (not math.isclose(row['unweighted_fixed_composite_loss'],combined,rel_tol=1e-6,abs_tol=1e-7) or
                not math.isclose(row['weighted_loss'],row['unweighted_fixed_composite_loss']*row['importance_weight'],
                                 rel_tol=1e-6,abs_tol=1e-7)):
            raise ValueError('Actual fixed composite/importance arithmetic differs')
    source=torch.load(D.ROOT/D.CONFIG['frozen_local_source'],map_location='cpu',weights_only=True)['model']
    model=D.F.A.new_model('cpu'); named=list(model.named_parameters()); final_exposure=None
    for epoch,step in enumerate((1536,3072,4608)):
        path=D.F.OUT/('step_'+str(step).zfill(6)+'.pt'); report=fitted['epochs'][epoch]
        if (report['epoch']!=epoch+1 or report['steps']!=1536 or report['counts']!=dict(expected_counts) or
                report['checkpoint']!=str(path.relative_to(D.ROOT)) or not path.stat().st_size):
            raise ValueError('Actual saved epoch/checkpoint exposure differs')
        state=torch.load(path,map_location='cpu',weights_only=True)
        if (state['arm']!='A40' or state['configuration']!=D.CONFIG or state['protocol']!=D.PROTOCOL or
                state['step']!=step or state['revision']!=D.F.A.R.C.I.F.REVISION or
                sum(v.numel() for v in state['model'].values())!=613348 or
                any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values()) or
                state['frozen_local_source_checkpoint']!=D.CONFIG['frozen_local_source'] or
                set(state['frozen_local_source_model'])!=set(source) or
                any(not torch.equal(v,source[k]) for k,v in state['frozen_local_source_model'].items())):
            raise ValueError('Finite fixed student and exact unchanged original A3 source required')
        model.load_state_dict(state['model'],strict=True)
        opt=state['optimizer']; groups=opt['param_groups']
        if (len(groups)!=1 or len(groups[0]['params'])!=len(named) or
                len(set(groups[0]['params']))!=len(named) or groups[0]['lr']!=D.CONFIG['lr'] or
                groups[0]['weight_decay']!=D.CONFIG['weight_decay'] or list(groups[0]['betas'])!=D.CONFIG['betas'] or
                groups[0]['eps']!=D.CONFIG['eps']):
            raise ValueError('Actual fixed student-only AdamW configuration differs')
        exposure=expected_optimizer_steps([k for k,_ in named],draws,records,epoch+1)
        expected_ids={ident for ident,(name,_) in zip(groups[0]['params'],named) if exposure[name]}
        if set(opt['state'])!=expected_ids:
            raise ValueError('Actual student branch optimizer population differs')
        for ident,(name,param) in zip(groups[0]['params'],named):
            if not exposure[name]: continue
            value=opt['state'][ident]
            if (float(value['step'])!=exposure[name] or
                    any(v.shape!=param.shape or v.dtype!=torch.float32 or not bool(torch.isfinite(v).all())
                        for v in (value['exp_avg'],value['exp_avg_sq']))):
                raise ValueError('Actual finite optimizer exposure differs: '+name)
        final_exposure=exposure; del state,opt
    if fitted['optimizer_steps_per_parameter']!=final_exposure:
        raise ValueError('Actual branch-step summary differs from checkpoint and original pool')
    return dict(expected_counts),final_exposure


@torch.inference_mode()
def main(fit_job,readout_job):
    began=time.perf_counter(); D.F.A.A.require_jobs([fit_job,readout_job])
    report=D.R.completed(D.OUT,readout_job)
    if (report['configuration']!=D.CONFIG or report['protocol']!=D.PROTOCOL or report['source_fit_job']!=fit_job or
            report['source_preflight_job']!=D.F.PREFLIGHT_JOB or report['windows']!=1491 or
            report['natural_candidates']!=16321 or report['reference_injected_predictions']!=0 or
            report['local_logit_receipts']!=1491 or report['optimizer_steps_this_run']!=0 or
            any(report[k] is not True for k in ('natural_score_file_closed_before_saved_label_arrays_open',
                'local_logit_files_closed_before_saved_label_arrays_open','original_natural_geometry_order_sources_unchanged',
                'original_saved_labels_reused','model_weights_unchanged','frozen_local_weights_unchanged',
                'all_parameters_frozen_and_gradients_absent')) or any(report[k] is not False for k in
                ('label_catalog_opened','new_candidate_search','new_backbone_forward','new_R4_labelling',
                 'DEV_used','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual complete unchanged label-free natural A40 readout required')
    replay,exposure=verify_fit(fit_job)
    _,free,previous,_,baseline=D.B.population()
    saved=D.R.C.R.rows(D.OUT/'natural_scores.jsonl'); labelled=D.R.C.R.rows(D.OUT/'windows.jsonl')
    if len(saved)!=1491 or len(labelled)!=1491:
        raise ValueError('All1491 original natural score/label records required')
    maximum=0.
    for ordinal,(row,pool,old) in enumerate(zip(saved,free,previous)):
        path=D.OUT/'local_logits'/('window_'+str(ordinal).zfill(6)+'.pt')
        if (row['window']!=pool['window'] or row['chains']!=pool['chains'] or
                row['decoder_sources']!=pool['decoder_sources'] or row['local_logit_artifact']!=str(path.relative_to(D.ROOT))):
            raise ValueError('Original saved natural identity/order/source differs')
        verify_components(row)
        receipt=torch.load(path,map_location='cpu',weights_only=True)
        logits=receipt['local_logits']
        if (receipt['window']!=row['window'] or receipt['source_checkpoint']!=D.CONFIG['frozen_local_source'] or
                receipt['revision']!=D.F.A.R.C.I.F.REVISION or list(logits.shape)!=[24576,15] or
                logits.dtype!=torch.float32 or not bool(torch.isfinite(logits).all())):
            raise ValueError('Actual original frozen oriented FP32 local input receipt differs')
        actual=local_path_scores(logits,row['chains'],row['window']['valid_bases'])
        expected=torch.tensor(old['variants']['A31_A3']['scores'],dtype=torch.float64)
        torch.testing.assert_close(actual,expected,rtol=1e-12,atol=1e-8)
        maximum=max(maximum,float((actual-expected).abs().max()) if len(actual) else 0.)
        del receipt,logits,actual,expected
    original_labels=D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    prior=D.R.C.R.rows(D.F.A.M.D.OUT/'natural_scores.jsonl')
    closed=D.R.completed(D.F.A.M.D.CHECK,D.F.A.M.I.CHECK_JOB)
    rows,counts,references=D.paired_readout(saved,previous,original_labels,prior)
    if (rows!=labelled or counts!=report['selected_counts_per_species'] or references!=report['eligible_reference_pairs'] or
            sum(v['eligible'] for v in references.values())!=2916 or sum(v['five_or_more'] for v in references.values())!=1086 or
            any(counts[s]['A31_A3']!=baseline['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A39']!=closed['selected_counts_per_species'][s]['A39'] for s in counts)):
        raise ValueError('Actual protected labels/DP/full denominators/stronger baseline and closed A39 differ')
    checks=D.gate(counts)
    if checks!=report['predeclared_TRAIN_gate'] or all(checks.values())!=report['TRAIN_preparation_supported']:
        raise ValueError('Actual original five-condition arithmetic differs')
    D.R.C.R.write_summary(D.CHECK,dict(completed=True,experiment='M28-A40-ROLE-EVIDENCE-CONTEXT-RESULT-CHECK',
        protocol=D.PROTOCOL,configuration=D.CONFIG,source_jobs=dict(fit=fit_job,readout=readout_job),
        windows=1491,natural_candidates=16321,reference_injected_predictions=0,
        actual_original_draw_weights_counts_and_4608_updates_verified=True,
        actual_saved_fixed_composite_loss_and_importance_arithmetic_verified=True,
        actual_three_student_checkpoints_and_branch_optimizer_exposure_verified=True,
        actual_three_frozen_source_states_equal_original_A3_and_excluded_from_optimizer_verified=True,
        actual_saved_frozen_A3_role_inputs_and_original_A31_scores_verified=True,
        actual_saved_natural_identity_and_FP32_score_components_verified=True,
        actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed=True,
        original_replayed_per_pass_counts=replay,actual_optimizer_steps_per_parameter=exposure,
        local_logit_receipts=1491,max_original_A31_local_score_error=maximum,
        selected_counts_per_species=counts,eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks,TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False,independent_neural_forward_replay=False,optimizer_steps_this_run=0,
        new_candidate_search=False,new_R4_labelling=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='Exact saved FP32 role input/component and original TRAIN decision verification; '
            'not independent neural replay, biological causality, genome gene F1 or independent generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-job',required=True); parser.add_argument('--readout-job',required=True)
    args=parser.parse_args(); torch.set_num_threads(2); main(args.fit_job,args.readout_job)
