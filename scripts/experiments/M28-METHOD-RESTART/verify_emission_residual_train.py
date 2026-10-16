#!/usr/bin/env python3
"""Actual A41 affine exposure, frozen states, saved potentials and TRAIN decisions."""
import argparse
from collections import Counter
import math
import time
import torch
import readout_emission_residual_train as D
from src.m28.local_path_score import local_path_scores


def verify_optimizer(opt,named,step):
    groups=opt['param_groups']
    if (len(groups)!=1 or len(groups[0]['params'])!=len(named) or
            len(set(groups[0]['params']))!=len(named) or groups[0]['lr']!=D.CONFIG['lr'] or
            groups[0]['weight_decay']!=D.CONFIG['weight_decay'] or list(groups[0]['betas'])!=D.CONFIG['betas'] or
            groups[0]['eps']!=D.CONFIG['eps'] or set(opt['state'])!=set(groups[0]['params'])):
        raise ValueError('Fixed residual-only AdamW population and configuration required')
    for ident,(name,param) in zip(groups[0]['params'],named):
        saved=opt['state'][ident]
        if (float(saved['step'])!=step or
                any(v.shape!=param.shape or v.dtype!=torch.float32 or not bool(torch.isfinite(v).all())
                    for v in (saved['exp_avg'],saved['exp_avg_sq']))):
            raise ValueError('Actual finite affine optimizer exposure differs: '+name)


@torch.inference_mode()
def verify_fit(fit_job):
    fitted=D.R.completed(D.F.OUT,fit_job); D.validate_fit(fitted)
    draws,_,_,records,_,_,_=D.F.A.A.population()
    training=D.R.C.R.rows(D.F.OUT/'training.jsonl')
    if len(training)!=4608 or len(fitted['epochs'])!=3:
        raise ValueError('Actual three complete original passes required')
    replay=Counter()
    for draw in draws: replay.update(records[draw['window_id']]['counts'])
    fields=('weighted_loss','unweighted_mean_known_BCE','unweighted_mean_known_set_NLL',
        'unweighted_fixed_composite_loss','log_partition','gold_gain','gradient_norm_before_clip')
    for i,row in enumerate(training):
        epoch,ordinal=divmod(i,1536); draw=draws[ordinal]; rec=records[draw['window_id']]
        if (row['step']!=i+1 or row['epoch']!=epoch+1 or row['draw_ordinal']!=ordinal or
                row['window_id']!=draw['window_id'] or row['species']!=draw['species'] or
                row['importance_weight']!=draw['importance_weight_uniform_within_stratum'] or
                row['frozen_source_gradients_absent'] is not True or any(row[k]!=v for k,v in rec['counts'].items()) or
                any(not math.isfinite(row[k]) for k in fields) or
                (rec['counts']['zero_known_windows'] and any(row[k]!=0 for k in fields))):
            raise ValueError('Actual original draw, weight, protected update or finite counts differ')
        if (not math.isclose(row['unweighted_fixed_composite_loss'],
                (row['unweighted_mean_known_BCE']+row['unweighted_mean_known_set_NLL'])/2,rel_tol=1e-12,abs_tol=1e-8) or
                not math.isclose(row['weighted_loss'],row['unweighted_fixed_composite_loss']*row['importance_weight'],
                                 rel_tol=1e-12,abs_tol=1e-8)):
            raise ValueError('Actual FP64 fixed composite and importance arithmetic differs')
    original=torch.load(D.ROOT/D.CONFIG['frozen_source'],map_location='cpu',weights_only=True)['model']
    model=D.F.A.new_model('cpu'); named=list(model.named_parameters())
    for epoch,step in enumerate((1536,3072,4608)):
        path=D.F.OUT/('step_'+str(step).zfill(6)+'.pt'); saved=fitted['epochs'][epoch]
        actual_mean=sum(r['weighted_loss'] for r in training[epoch*1536:(epoch+1)*1536])/1536
        if (saved['epoch']!=epoch+1 or saved['steps']!=1536 or saved['counts']!=dict(replay) or
                saved['checkpoint']!=str(path.relative_to(D.ROOT)) or not path.stat().st_size or
                not math.isclose(saved['mean_weighted_loss'],actual_mean,rel_tol=1e-12,abs_tol=1e-8)):
            raise ValueError('Actual saved epoch and checkpoint exposure differs')
        state=torch.load(path,map_location='cpu',weights_only=True)
        if (state['arm']!='A41' or state['configuration']!=D.CONFIG or state['protocol']!=D.PROTOCOL or
                state['step']!=step or state['revision']!=D.R.C.I.F.REVISION or
                sum(v.numel() for v in state['model'].values())!=1806 or
                any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values()) or
                state['frozen_source_checkpoint']!=D.CONFIG['frozen_source'] or
                set(state['frozen_source_model'])!=set(original) or
                any(not torch.equal(v,original[k]) for k,v in state['frozen_source_model'].items())):
            raise ValueError('Finite affine state and unchanged original frozen A3 source required')
        model.load_state_dict(state['model'],strict=True); verify_optimizer(state['optimizer'],named,step)
        del state
    exposure={k:4608 for k,_ in named}
    if fitted['optimizer_steps_per_parameter']!=exposure:
        raise ValueError('Actual two-parameter step summary differs')
    return dict(replay),exposure


def verify_potentials(row,receipt,original_scores):
    length=row['window']['valid_bases']; frozen=receipt['frozen_local_logits']; residual=receipt['residual_local_logits']
    if (receipt['window']!=row['window'] or receipt['source_checkpoint']!=D.CONFIG['frozen_source'] or
            receipt['revision']!=D.R.C.I.F.REVISION or row['null']!=0. or
            any(list(v.shape)!=[length,15] or v.dtype!=torch.float32 or not bool(torch.isfinite(v).all())
                for v in (frozen,residual)) or not bool((residual[:,0]==0).all())):
        raise ValueError('Actual oriented FP32 frozen and residual local receipts differ')
    base=local_path_scores(frozen,row['chains'],length); correction=local_path_scores(residual,row['chains'],length)
    expected=torch.tensor(original_scores,dtype=torch.float64)
    arrays={k:torch.tensor(row[k],dtype=torch.float64) for k in ('scores','frozen_path_gain','residual_path_gain')}
    if len(expected)!=len(row['chains']) or any(len(v)!=len(row['chains']) or not bool(torch.isfinite(v).all())
                                               for v in arrays.values()):
        raise ValueError('One finite FP64 potential per original candidate required')
    torch.testing.assert_close(base,expected,rtol=1e-12,atol=1e-8)
    torch.testing.assert_close(base,arrays['frozen_path_gain'],rtol=1e-12,atol=1e-8)
    torch.testing.assert_close(correction,arrays['residual_path_gain'],rtol=1e-12,atol=1e-8)
    # JSON preserves FP64 components; never round the final sum to FP32.
    torch.testing.assert_close(arrays['scores'],arrays['frozen_path_gain']+arrays['residual_path_gain'],rtol=0,atol=0)
    return float((base-expected).abs().max()) if len(base) else 0.


@torch.inference_mode()
def main(fit_job,readout_job):
    began=time.perf_counter(); D.F.A.A.require_jobs([fit_job,readout_job])
    report=D.R.completed(D.OUT,readout_job)
    if (report['configuration']!=D.CONFIG or report['protocol']!=D.PROTOCOL or report['source_fit_job']!=fit_job or
            report['source_preflight_job']!=D.F.PREFLIGHT_JOB or report['windows']!=1491 or
            report['natural_candidates']!=16321 or report['potential_receipts']!=1491 or
            report['reference_injected_predictions']!=0 or report['optimizer_steps_this_run']!=0 or
            any(report[k] is not True for k in ('all_scores_and_potential_files_closed_before_saved_labels_open',
                'original_natural_geometry_order_sources_unchanged','original_saved_labels_reused',
                'model_weights_unchanged','frozen_source_weights_unchanged','all_parameters_frozen_and_gradients_absent')) or
            any(report[k] is not False for k in ('label_catalog_opened','new_candidate_search','new_backbone_forward',
                'new_R4_labelling','DEV_used','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual complete unchanged label-free A41 readout required')
    replay,exposure=verify_fit(fit_job)
    _,free,previous,_,baseline=D.B.population()
    saved=D.R.C.R.rows(D.OUT/'natural_scores.jsonl'); labelled=D.R.C.R.rows(D.OUT/'windows.jsonl')
    if len(saved)!=1491 or len(labelled)!=1491:
        raise ValueError('All1491 original natural records required')
    maximum=0.
    for ordinal,(row,pool,old) in enumerate(zip(saved,free,previous)):
        path=D.OUT/'potentials'/('window_'+str(ordinal).zfill(6)+'.pt')
        if (row['window']!=pool['window'] or row['chains']!=pool['chains'] or
                row['decoder_sources']!=pool['decoder_sources'] or row['potential_artifact']!=str(path.relative_to(D.ROOT))):
            raise ValueError('Original natural identity, order and provenance differ')
        receipt=torch.load(path,map_location='cpu',weights_only=True)
        maximum=max(maximum,verify_potentials(row,receipt,old['variants']['A31_A3']['scores']))
        del receipt
    source_labels=D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    rows,counts,refs=D.paired_readout(saved,previous,source_labels)
    if (rows!=labelled or counts!=report['selected_counts_per_species'] or refs!=report['eligible_reference_pairs'] or
            sum(v['eligible'] for v in refs.values())!=2916 or sum(v['five_or_more'] for v in refs.values())!=1086 or
            any(counts[s]['A31_A3']!=baseline['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Actual protected labels, full DP counts and denominators differ')
    checks=D.gate(counts)
    if checks!=report['predeclared_TRAIN_gate'] or all(checks.values())!=report['TRAIN_preparation_supported']:
        raise ValueError('Actual original five-condition arithmetic differs')
    D.R.C.R.write_summary(D.CHECK,dict(completed=True,protocol=D.PROTOCOL,configuration=D.CONFIG,
        source_jobs=dict(fit=fit_job,readout=readout_job),windows=1491,natural_candidates=16321,potential_receipts=1491,
        actual_original_draw_weights_counts_and_4608_updates_verified=True,
        actual_saved_fixed_composite_loss_and_importance_arithmetic_verified=True,
        actual_three_affine_checkpoints_and_optimizer_exposure_verified=True,
        actual_three_frozen_source_states_equal_original_A3_and_excluded_from_optimizer_verified=True,
        actual_saved_FP32_local_potentials_and_FP64_gain_arithmetic_verified=True,
        actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed=True,
        original_replayed_per_pass_counts=replay,actual_optimizer_steps_per_parameter=exposure,
        max_original_A31_score_error=maximum,selected_counts_per_species=counts,eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=checks,TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False,independent_neural_forward_replay=False,optimizer_steps_this_run=0,
        new_candidate_search=False,new_R4_labelling=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='Saved source and residual potentials, optimizer exposure and TRAIN decisions; '
            'not independent replay of affine outputs from context, raw-FASTA gene F1 or generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-job',required=True); parser.add_argument('--readout-job',required=True)
    args=parser.parse_args(); torch.set_num_threads(2); main(args.fit_job,args.readout_job)
