#!/usr/bin/env python3
"""Native weighted disposable updates and strict raw-adapter routing for the frozen TRAIN proposal."""
import copy
import json
import os
import time
import torch
import fit_full_grammar_objective_contrast as C
import prepare_full_grammar_complete_core_draws as P
import check_full_grammar_complete_core_draws_support as S
import annotate_full_grammar_complete_core_draws as A

ROOT = P.ROOT
PREFLIGHT = ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-PREFLIGHT'
CONFIG = dict(C.F.CONFIG,**A.IDENTITY,maximum_new_fit_bytes=96*1024*1024)


def plan_and_population():
    S.population()
    supported = json.loads((S.AGGREGATE/'summary.json').read_text())
    if ((S.AGGREGATE/'STATUS').read_text().strip() != 'COMPLETED' or not supported['completed']
            or supported['job_id'] != A.IDENTITY['support_job'] or supported['configuration'] != S.CONFIG
            or supported['unique_events'] != 4107 or supported['prospective_draws'] != 4608
            or not supported['all_fixed_draws_retained'] or not supported['all_selected_exact_fine_events_feasible']
            or supported['optimizer_steps_this_run'] != 0 or supported['DEV_used'] or supported['sealed_data_used']):
        raise ValueError('All fixed complete-core-proposal events must be certified feasible')
    plan = list(P.P.rows(P.OUT/'training_plan.jsonl'))
    windows,strata,catalogs,original = P.P.load_population()
    eligibility = {r['window_id']:r for r in P.P.rows(P.C.OUT/'window_eligibility.jsonl')}
    probabilities = {key:P.proposal(pool,eligibility) for key,pool in strata.items()}
    replay = P.schedule(original,strata,probabilities)
    if len(plan) != len(replay) or any({k:row[k] for k in expected} != expected for row,expected in zip(plan,replay)):
        raise ValueError('Fixed complete-core seed1 draw schedule or proposal probabilities changed')
    for row in plan:
        denominator = max(1,row['fine_positions']+row['coarse_positions'])
        iw = 1/(row['stratum_size']*row['stratum_probability'])
        if (row['importance_weight'] != iw or row['known_observation_denominator'] != denominator
                or row['objective_multiplier'] != iw/denominator):
            raise ValueError('Proposal importance or known-row denominator differs')
    return plan,windows,catalogs


def identity_fixture():
    state = dict(configuration=CONFIG,step=9216,epoch=6,protocol=A.PROTOCOL,encoder_revision=C.Q.F.REVISION)
    A.validate_checkpoint(state)
    for key,value in (('uniform_mixture_probability',1.),('unique_windows',4219),('epsilon',.2)):
        altered = copy.deepcopy(state)
        altered['configuration'][key] = value
        try:
            A.validate_checkpoint(altered)
        except ValueError:
            pass
        else:
            raise AssertionError('Raw adapter accepted changed '+key)
    return 'complete_core_identity_accepts_fixed_and_rejects_mixture_population_or_likelihood_drift'


def main():
    began = time.perf_counter()
    fixture = identity_fixture()
    plan,windows,catalogs = plan_and_population()
    selected = [next(i for i,r in enumerate(plan[:1536]) if r['species'] == species) for species in C.Q.ALLOW]
    C.write(PREFLIGHT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=A.PROTOCOL,selected_additional_steps=[i+1 for i in selected],support_job=A.IDENTITY['support_job']))
    state = C.starting_state()
    sequences,descriptors = C.P.read_train_DNA()
    backbone,tokenizer = C.G.load_native()
    records = []
    for ordinal in selected:
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome,geometry,context,offset,targets = C.example(w,expected,sequences,catalogs)
        model,optimizer = C.load_model_optimizer(state)
        joined = C.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
        raw,_ = C.objective('STRUCTURED_CONTINUE',context,joined,offset,w['valid_bases'],targets)
        total,norm = C.update(model,optimizer,backbone,raw,expected['objective_multiplier'])
        if (any(int(s['step']) != 4609 for s in optimizer.state.values())
                or all(torch.equal(p.cpu(),state['model'][k]) for k,p in model.state_dict().items())):
            raise ValueError('Independent disposable original-AdamW update4609 required')
        fixture_state = dict(configuration=CONFIG,protocol=A.PROTOCOL,step=9216,epoch=6,
            encoder_revision=C.Q.F.REVISION,model=model.state_dict())
        inference_model = A.model_from_state(fixture_state,'cuda')
        model.eval()
        with torch.inference_mode():
            trained = C.potentials(model,backbone,tokenizer,chromosome,geometry,w['strand'])
            restored = C.potentials(inference_model,backbone,tokenizer,chromosome,geometry,w['strand'])
        if not torch.equal(trained,restored):
            raise ValueError('Actual native logits differ after strict complete-core raw-adapter load')
        row = dict(additional_step=ordinal+1,window_id=w['id'],species=w['species'],
            raw_NLL=float(raw.detach()),raw_NLL_dtype=str(raw.dtype),weighted_total=total,
            importance_weight=expected['importance_weight'],objective_multiplier=expected['objective_multiplier'],
            known_observation_denominator=expected['known_observation_denominator'],gradient_norm_before_clip=norm,
            initial_optimizer_step=4608,final_disposable_optimizer_step=4609,
            actual_native_logits_raw_adapter_exact=True)
        records.append(row)
        print(json.dumps(row,allow_nan=False),flush=True)
        del model,optimizer,inference_model,joined,raw,trained,restored,fixture_state
    C.write(PREFLIGHT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG,protocol=A.PROTOCOL,support_job=A.IDENTITY['support_job'],fixture=fixture,
        records=records,original_TRAIN_DNA_sources=descriptors,disposable_updates_only=True,
        optimizer_steps_this_run=2,formal_training_started=False,checkpoints_written=False,
        raw_adapter_fixture_metadata_not_formal9216=True,actual_native_logits_raw_adapter_exact=True,
        saved_source_weights_and_optimizer_unchanged=True,pretrained_encoder_frozen_no_grad=True,
        cached_features_used=False,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
