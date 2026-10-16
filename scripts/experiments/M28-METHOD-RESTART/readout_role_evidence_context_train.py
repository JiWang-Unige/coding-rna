#!/usr/bin/env python3
"""Complete label-free A40 natural scores, then original saved labels join."""
import argparse
from collections import Counter
import json
import time
import torch
import fit_role_evidence_context as F

ROOT,CONFIG,PROTOCOL=F.ROOT,F.CONFIG,F.PROTOCOL
R,B=F.A.M.R,F.A.A.B
OUT=ROOT/'outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-TRAIN-READOUT'
CHECK=ROOT/'outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-RESULT-CHECK'
COMPONENTS=('additive','nonadditive','exon_additive','junction_additive')


def validate_fit(value,state=False):
    if (value['arm']!='A40' or value['configuration']!=CONFIG or value['protocol']!=PROTOCOL or
            value['step' if state else 'optimizer_steps']!=4608 or
            value['revision' if state else 'feature_revision']!=F.A.R.C.I.F.REVISION):
        raise ValueError('Only fixed complete A40 last004608 identity allowed')
    if not state and (value['completed'] is not True or value['preflight_job']!=F.PREFLIGHT_JOB or
            value['trainable_parameters']!=613348 or value['initialization']!=CONFIG['initialization'] or
            value['final_checkpoint']!=str((F.OUT/'step_004608.pt').relative_to(ROOT)) or
            value['frozen_local_source_checkpoint']!=CONFIG['frozen_local_source'] or
            value['frozen_local_source_optimizer_steps']!=0 or value['reference_injected_training_candidates']!=0 or
            any(value[k] is not True for k in ('frozen_supervised_A3_checkpoint_read',
                'frozen_local_source_weights_unchanged','frozen_local_source_gradients_absent')) or
            any(value[k] is not False for k in ('student_supervised_checkpoint_read','generator_decoding',
                'generator_updated','new_candidate_search','new_backbone_forward','new_R4_labelling',
                'DEV_used_for_selection','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual unchanged complete natural-only A40 fit required')


def gate(counts):
    checks={}
    for species,values in counts.items():
        old,new=values['A31_A3'],values['A40']
        checks[species+'_natural_selected_correct_strictly_increased']=new['selected_correct']>old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased']=new['selected_known_nonmatch']<=old['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased']=(
        sum(v['A40']['selected_five_or_more_correct'] for v in counts.values())>
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def paired_readout(saved,previous,labelled,prior):
    if not len(saved)==len(previous)==len(labelled)==len(prior):
        raise ValueError('Complete original four-way ordered windows required')
    totals={s:{a:Counter() for a in ('A31_A3','A39','A40')} for s in R.C.R.ALLOW}
    refs={s:Counter(eligible=0,five_or_more=0) for s in totals}; rows=[]
    for new,old,lab,last in zip(saved,previous,labelled,prior):
        baseline=old['variants']['A31_A3']
        if (any(v['window']!=new['window'] for v in (old,lab,last)) or
                any(v['reference_used_in_scoring'] is not False for v in (new,old,last)) or
                new['window']['split']!='train' or new['chains']!=baseline['chains'] or new['chains']!=last['chains']):
            raise ValueError('Unchanged ordered reference-free same-natural TRAIN identity required')
        labels=lab['labels']['A31_A3']
        if any(y not in (-1,0,1) for y in labels): raise ValueError('Original protected labels required')
        counts=dict(A31_A3=R.selected_counts(baseline,labels),A39=R.selected_counts(last,labels),
                    A40=R.selected_counts(new,labels))
        if counts['A31_A3']!=lab['counts']['A31_A3']:
            raise ValueError('Original stronger-baseline selection differs')
        species=new['window']['species']
        for arm,value in counts.items(): totals[species][arm].update(value)
        refs[species].update(eligible=lab['eligible_references'],five_or_more=lab['eligible_five_or_more'])
        rows.append(dict(window=new['window'],labels=labels,counts=counts,
            eligible_references=lab['eligible_references'],eligible_five_or_more=lab['eligible_five_or_more']))
    return rows,{s:{a:dict(c) for a,c in v.items()} for s,v in totals.items()},{s:dict(v) for s,v in refs.items()}


@torch.inference_mode()
def main(fit_job):
    F.A.A.require_jobs([fit_job,F.PREFLIGHT_JOB]); F.ready()
    fitted=R.completed(F.OUT,fit_job); validate_fit(fitted)
    state=torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True); validate_fit(state,True)
    student=F.A.new_model('cuda').eval().requires_grad_(False); student.load_state_dict(state['model'],strict=True)
    source=F.A.local_source()
    if (state['frozen_local_source_checkpoint']!=CONFIG['frozen_local_source'] or
            set(state['frozen_local_source_model'])!=set(source.state_dict()) or
            any(not torch.equal(v.cpu(),state['frozen_local_source_model'][k]) for k,v in source.state_dict().items())):
        raise ValueError('Readout frozen local source differs from actual fit')
    del state
    before={name:{k:v.detach().clone() for k,v in model.state_dict().items()}
            for name,model in (('student',student),('source',source))}
    draws,free,previous,index,baseline=B.population()
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); (OUT/'local_logits').mkdir()
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal,pool in enumerate(free):
            w=pool['window']; cached=torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            logits=F.A.local_logits(source,features,dna,w['valid_bases'])
            score=student(features,dna,[w['valid_bases']],pool['chains'],logits)
            if any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
                raise ValueError('Finite unchanged FP32 role-conditioned scores required')
            path=OUT/'local_logits'/('window_'+str(ordinal).zfill(6)+'.pt')
            with path.open('xb') as file:
                torch.save(dict(window=w,source_checkpoint=CONFIG['frozen_local_source'],
                    revision=F.A.R.C.I.F.REVISION,local_logits=logits.cpu()),file)
            values=score['logits'].cpu().tolist()
            row=dict(window=w,chains=pool['chains'],decoder_sources=pool['decoder_sources'],scores=values,
                selected=R.select_nonoverlapping(pool['chains'],values),
                **{k:score[k].cpu().tolist() for k in COMPONENTS},null=float(score['null']),
                local_logit_artifact=str(path.relative_to(ROOT)),reference_used_in_scoring=False)
            handle.write(json.dumps(row,allow_nan=False)+'\n')
            if (ordinal+1)%64==0 or ordinal+1==len(free):
                handle.flush(); print(json.dumps(dict(natural_windows=ordinal+1,total=1491)),flush=True)
            del cached,features,dna,logits,score
    for name,model in (('student',student),('source',source)):
        if (any(not torch.equal(v,before[name][k]) for k,v in model.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Readout changed model/source or received gradients')
    # No label array joins until every new natural score and local receipt is closed.
    saved=R.C.R.rows(OUT/'natural_scores.jsonl'); labelled=R.C.R.rows(R.READOUT/'windows.jsonl')
    prior=R.C.R.rows(F.A.M.D.OUT/'natural_scores.jsonl')
    old=R.completed(F.A.M.D.CHECK,F.A.M.I.CHECK_JOB)
    rows,counts,references=paired_readout(saved,previous,labelled,prior)
    if (len(rows)!=1491 or sum(len(v['chains']) for v in saved)!=16321 or
            sum(v['eligible'] for v in references.values())!=2916 or sum(v['five_or_more'] for v in references.values())!=1086 or
            any(counts[s]['A31_A3']!=baseline['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A39']!=old['selected_counts_per_species'][s]['A39'] for s in counts)):
        raise ValueError('Original full counts/denominators/stronger baseline and A39 control differ')
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in rows: handle.write(json.dumps(row,allow_nan=False)+'\n')
    checks=gate(counts)
    R.C.R.write_summary(OUT,dict(completed=True,experiment='M28-A40-ROLE-EVIDENCE-CONTEXT-TRAIN-READOUT',
        protocol=PROTOCOL,fit_protocol=PROTOCOL,configuration=CONFIG,source_fit_job=fit_job,
        source_preflight_job=F.PREFLIGHT_JOB,source_baseline_jobs=dict(readout='13398521',CPU='13398524'),
        source_previous_arm_jobs=dict(fit='13415276',readout='13415309',CPU='13415319'),
        checkpoint=fitted['final_checkpoint'],frozen_local_source_checkpoint=CONFIG['frozen_local_source'],
        windows=1491,original_draws=1536,natural_candidates=16321,reference_injected_predictions=0,
        local_logit_receipts=1491,natural_score_file_closed_before_saved_label_arrays_open=True,
        local_logit_files_closed_before_saved_label_arrays_open=True,original_natural_geometry_order_sources_unchanged=True,
        original_saved_labels_reused=True,label_catalog_opened=False,selected_counts_per_species=counts,
        eligible_reference_pairs=references,predeclared_TRAIN_gate=checks,TRAIN_preparation_supported=all(checks.values()),
        model_weights_unchanged=True,frozen_local_weights_unchanged=True,all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0,new_candidate_search=False,new_backbone_forward=False,new_R4_labelling=False,
        DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Complete in-sample original natural window/chain decisions with frozen supervised A3 role inputs; '
            'not raw-FASTA gene F1, matched-capacity causality or independent generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--fit-job',required=True)
    args=parser.parse_args(); torch.set_num_threads(2); main(args.fit_job)
