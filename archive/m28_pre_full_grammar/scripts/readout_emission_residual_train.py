#!/usr/bin/env python3
"""All label-free A41 potentials close before original saved labels join."""
import argparse
from collections import Counter
import json
import time
import torch
import fit_emission_residual as F

ROOT,R,CONFIG,PROTOCOL=F.ROOT,F.R,F.CONFIG,F.PROTOCOL
B=F.A.A.B
OUT=ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-TRAIN-READOUT'
CHECK=ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-RESULT-CHECK'
CHECKS=ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-READOUT-CHECKS'


def validate_fit(value,state=False):
    if (value['arm']!='A41' or value['configuration']!=CONFIG or value['protocol']!=PROTOCOL or
            value['step' if state else 'optimizer_steps']!=4608 or
            value['revision' if state else 'feature_revision']!=R.C.I.F.REVISION):
        raise ValueError('Fixed complete A41 last004608 required')
    if not state and (value['completed'] is not True or value['preflight_job']!=F.PREFLIGHT_JOB or
            value['trainable_parameters']!=1806 or value['initialization']!=CONFIG['initialization'] or
            value['final_checkpoint']!=str((F.OUT/'step_004608.pt').relative_to(ROOT)) or
            value['frozen_source_checkpoint']!=CONFIG['frozen_source'] or value['frozen_source_optimizer_steps']!=0 or
            value['reference_injected_training_candidates']!=0 or
            any(value[k] is not True for k in ('frozen_source_weights_unchanged','frozen_source_gradients_absent')) or
            any(value[k] is not False for k in ('generator_decoding','generator_updated','new_candidate_search',
                'new_backbone_forward','new_R4_labelling','DEV_used_for_selection','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual unchanged natural-only A41 fit required')


def gate(counts):
    checks={}
    for species,arms in counts.items():
        new,old=arms['A41'],arms['A31_A3']
        checks[species+'_natural_selected_correct_strictly_increased']=new['selected_correct']>old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased']=new['selected_known_nonmatch']<=old['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased']=(
        sum(v['A41']['selected_five_or_more_correct'] for v in counts.values())>
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def paired_readout(saved,previous,labelled):
    if not len(saved)==len(previous)==len(labelled):
        raise ValueError('Complete original ordered three-way windows required')
    counts={s:{a:Counter() for a in ('A31_A3','A41')} for s in R.C.R.ALLOW}
    refs={s:Counter(eligible=0,five_or_more=0) for s in counts}; rows=[]
    for new,old,lab in zip(saved,previous,labelled):
        baseline=old['variants']['A31_A3']
        if (new['window']!=old['window'] or new['window']!=lab['window'] or
                new['window']['split']!='train' or new['chains']!=baseline['chains'] or
                new['reference_used_in_scoring'] is not False or old['reference_used_in_scoring'] is not False):
            raise ValueError('Original ordered reference-free natural identity differs')
        labels=lab['labels']['A31_A3']
        if any(y not in (-1,0,1) for y in labels): raise ValueError('Original protected labels required')
        local=dict(A31_A3=R.selected_counts(baseline,labels),A41=R.selected_counts(new,labels))
        if local['A31_A3']!=lab['counts']['A31_A3']:
            raise ValueError('Original stronger-baseline counts differ')
        species=new['window']['species']
        for arm,value in local.items(): counts[species][arm].update(value)
        refs[species].update(eligible=lab['eligible_references'],five_or_more=lab['eligible_five_or_more'])
        rows.append(dict(window=new['window'],labels=labels,counts=local,
            eligible_references=lab['eligible_references'],eligible_five_or_more=lab['eligible_five_or_more']))
    return rows,{s:{a:dict(c) for a,c in v.items()} for s,v in counts.items()},{s:dict(v) for s,v in refs.items()}


def prepare():
    _,free,previous,_,baseline=B.population()
    labelled=R.C.R.rows(R.READOUT/'windows.jsonl')
    placeholder=[dict(window=x['window'],reference_used_in_scoring=False,**x['variants']['A31_A3']) for x in previous]
    rows,counts,refs=paired_readout(placeholder,previous,labelled)
    if (len(rows)!=1491 or sum(len(x['chains']) for x in free)!=16321 or
            sum(x['eligible'] for x in refs.values())!=2916 or sum(x['five_or_more'] for x in refs.values())!=1086 or
            any(counts[s]['A31_A3']!=baseline['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A41']!=counts[s]['A31_A3'] for s in counts)):
        raise ValueError('Original full population, label join and denominators differ')
    R.C.R.write_summary(CHECKS,dict(completed=True,protocol=PROTOCOL,configuration=CONFIG,
        windows=1491,natural_candidates=16321,eligible_reference_pairs=refs,
        original_paired_population_ready=True,placeholder_is_original_baseline_not_A41_result=True,
        new_fitted_checkpoint_opened=False,feature_tensors_opened=False,model_or_backbone_forward=False,
        optimizer_steps_this_run=0,new_R4_labelling=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False))


@torch.inference_mode()
def main(fit_job,checks_job):
    F.A.A.require_jobs([fit_job,checks_job]); F.ready()
    prepared=R.completed(CHECKS,checks_job)
    if (prepared['configuration']!=CONFIG or prepared['protocol']!=PROTOCOL or
            prepared['original_paired_population_ready'] is not True):
        raise ValueError('Actual original paired population readiness required')
    fitted=R.completed(F.OUT,fit_job); validate_fit(fitted)
    state=torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True); validate_fit(state,True)
    model=F.A.new_model('cuda').eval().requires_grad_(False); model.load_state_dict(state['model'],strict=True)
    source=F.source_model()
    if (state['frozen_source_checkpoint']!=CONFIG['frozen_source'] or
            set(state['frozen_source_model'])!=set(source.state_dict()) or
            any(not torch.equal(v.cpu(),state['frozen_source_model'][k]) for k,v in source.state_dict().items())):
        raise ValueError('Actual fit frozen source differs from original A3')
    del state
    before={name:{k:v.detach().clone() for k,v in obj.state_dict().items()}
            for name,obj in (('residual',model),('source',source))}
    draws,free,previous,index,baseline=B.population()
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); (OUT/'potentials').mkdir()
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal,pool in enumerate(free):
            w=pool['window']; length=w['valid_bases']
            cached=torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            hidden,local=F.local_evidence(source,features,dna,length)
            score=model(hidden,local,pool['chains'],length)
            if any(not bool(torch.isfinite(v).all()) for v in score.values()):
                raise ValueError('Finite unchanged frozen and learned potentials required')
            path=OUT/'potentials'/('window_'+str(ordinal).zfill(6)+'.pt')
            with path.open('xb') as file:
                torch.save(dict(window=w,source_checkpoint=CONFIG['frozen_source'],revision=R.C.I.F.REVISION,
                    frozen_local_logits=local[:length].cpu(),residual_local_logits=score['residual_local_logits'].cpu()),file)
            values=score['logits'].cpu().tolist()
            handle.write(json.dumps(dict(window=w,chains=pool['chains'],decoder_sources=pool['decoder_sources'],
                scores=values,selected=R.select_nonoverlapping(pool['chains'],values),
                frozen_path_gain=score['frozen_path_gain'].cpu().tolist(),
                residual_path_gain=score['residual_path_gain'].cpu().tolist(),null=float(score['null']),
                potential_artifact=str(path.relative_to(ROOT)),reference_used_in_scoring=False),allow_nan=False)+'\n')
            if (ordinal+1)%64==0 or ordinal+1==len(free):
                handle.flush(); print(json.dumps(dict(natural_windows=ordinal+1,total=1491)),flush=True)
            del cached,features,dna,hidden,local,score
    for name,obj in (('residual',model),('source',source)):
        if (any(not torch.equal(v,before[name][k]) for k,v in obj.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in obj.parameters())):
            raise ValueError('Readout modified weights or received gradients')
    saved=R.C.R.rows(OUT/'natural_scores.jsonl'); labelled=R.C.R.rows(R.READOUT/'windows.jsonl')
    rows,counts,refs=paired_readout(saved,previous,labelled)
    if (len(rows)!=1491 or sum(len(x['chains']) for x in saved)!=16321 or
            sum(x['eligible'] for x in refs.values())!=2916 or sum(x['five_or_more'] for x in refs.values())!=1086 or
            any(counts[s]['A31_A3']!=baseline['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Complete original stronger baseline and denominators required')
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in rows: handle.write(json.dumps(row,allow_nan=False)+'\n')
    checks=gate(counts)
    R.C.R.write_summary(OUT,dict(completed=True,protocol=PROTOCOL,configuration=CONFIG,
        source_fit_job=fit_job,source_checks_job=checks_job,source_preflight_job=F.PREFLIGHT_JOB,
        checkpoint=fitted['final_checkpoint'],windows=1491,original_draws=len(draws),natural_candidates=16321,
        reference_injected_predictions=0,potential_receipts=1491,selected_counts_per_species=counts,
        eligible_reference_pairs=refs,predeclared_TRAIN_gate=checks,TRAIN_preparation_supported=all(checks.values()),
        all_scores_and_potential_files_closed_before_saved_labels_open=True,
        original_natural_geometry_order_sources_unchanged=True,original_saved_labels_reused=True,
        model_weights_unchanged=True,frozen_source_weights_unchanged=True,all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0,label_catalog_opened=False,new_candidate_search=False,new_backbone_forward=False,
        new_R4_labelling=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Complete fixed-pool TRAIN potentials, not raw-FASTA gene F1 or independent generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',action='store_true'); parser.add_argument('--fit-job'); parser.add_argument('--checks-job')
    args=parser.parse_args(); torch.set_num_threads(2)
    if args.prepare: prepare()
    elif args.fit_job and args.checks_job: main(args.fit_job,args.checks_job)
    else: parser.error('Actual fit and preparation job IDs required')
