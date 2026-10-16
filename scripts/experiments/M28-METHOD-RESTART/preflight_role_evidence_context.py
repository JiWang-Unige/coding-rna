#!/usr/bin/env python3
"""A40 zero-update candidate-role tokens on unchanged natural TRAIN geometry."""
import math
import time
import torch
import preflight_context_chain as A
import diagnose_null_anchored_set_attrition as M
from src.m28.role_evidence_context import RoleEvidenceContext
from src.m28.local_path_score import local_path_scores
from src.m28.null_anchored_chain_set import null_anchored_chain_set_loss

ROOT, R = A.ROOT, A.R
OUT = ROOT/'outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a40_role_evidence_context_protocol.md'
CONFIG = {**M.D.CONFIG, 'model':'fresh_joint_context_with_candidate_role_conditioned_exon_junction_tokens',
    'trainable_parameters':613348, 'head_parameters':249092,
    'initialization':'one_seed0_stream_SharedContext_then_expanded_role_conditioned_head_no_old_student_weights',
    'frozen_local_source':'outputs/M28-A3-FIRST-FIT/step_004608.pt',
    'frozen_local_source_parameters':649457,
    'role_evidence':'FP32_candidate_role_minus_intergenic_margin_3_exon_5_junction_token_columns',
    'local_source_updated':False}


def new_model(device):
    torch.manual_seed(CONFIG['seed'])
    return RoleEvidenceContext(CONFIG['feature_dim'],CONFIG['hidden']).to(device)


def local_source():
    model = R.C.load_generator().to('cuda').eval().requires_grad_(False)
    if sum(p.numel() for p in model.parameters()) != CONFIG['frozen_local_source_parameters']:
        raise ValueError('Original frozen A3 local source structure differs')
    return model


@torch.no_grad()
def local_logits(source,features,dna,length):
    hidden,_ = source.shared(features,dna,[length])
    logits = source.segmentation(hidden[0])
    if logits.dtype != torch.float32 or not bool(torch.isfinite(logits[:length]).all()):
        raise ValueError('Original frozen FP32 local evidence required')
    return logits


def closed_a39():
    A.require_jobs([M.I.FIT_JOB,M.I.READOUT_JOB,M.I.CHECK_JOB,'13417784'])
    fitted=M.R.completed(M.D.F.OUT,M.I.FIT_JOB)
    report=M.R.completed(M.D.OUT,M.I.READOUT_JOB)
    checked=M.R.completed(M.D.CHECK,M.I.CHECK_JOB)
    gates=M.validate_completed(fitted,report,checked)
    attrition=M.R.completed(M.OUT,'13417784')
    if (all(gates.values()) or attrition['scientific_status']!='A39_TRAIN_FAIL_closed_no_rescue' or
            attrition['selected_counts_per_species']!=report['selected_counts_per_species'] or
            attrition['eligible_reference_pairs']!=report['eligible_reference_pairs']):
        raise ValueError('Actual A39 closed FAIL must remain separate from the new representation')


def main():
    closed_a39()
    draws,index,windows,records,source_records,baseline,references=A.population()
    student=new_model('cuda').train(); source=local_source()
    if (sum(p.numel() for p in student.parameters())!=CONFIG['trainable_parameters'] or
            sum(p.numel() for p in student.shared.parameters())!=CONFIG['context_parameters'] or
            sum(p.numel() for p in student.head.parameters())!=CONFIG['head_parameters']):
        raise ValueError('Declared fresh role-conditioned structure differs')
    states={name:{k:v.detach().clone() for k,v in model.state_dict().items()}
            for name,model in (('student',student),('source',source))}
    old={row['window']['id']:row['variants']['A31_A3']
         for row in R.C.R.rows(M.R.READOUT/'natural_scores.jsonl')}
    weights={d['window_id']:d['importance_weight_uniform_within_stratum'] for d in draws}
    began=time.perf_counter(); torch.cuda.reset_peak_memory_stats(); receipts=[]
    for witness in R.P.select_gradient_windows(records):
        ident=witness['window_id']; w=windows[ident]; rec=records[ident]
        chains,labels=[r['chain'] for r in rec['candidates']],[r['label'] for r in rec['candidates']]
        applicable=any(y>=0 and len(c)>1 for c,y in zip(chains,labels))
        started=time.perf_counter()
        cached=torch.load(ROOT/index[ident]['artifact'],map_location='cpu',weights_only=True)
        features,dna=R.C.R.Q.window_inputs(cached,w,'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Original encoder inputs must remain detached')
        logits=local_logits(source,features,dna,w['valid_bases'])
        reconstructed=local_path_scores(logits,chains,w['valid_bases'])
        expected=torch.tensor(old[ident]['scores'],dtype=torch.float64,device='cuda')
        torch.testing.assert_close(reconstructed,expected,rtol=1e-12,atol=1e-8)
        student.zero_grad(set_to_none=True)
        score=student(features,dna,[w['valid_bases']],chains,logits)
        if any(v.dtype!=torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
            raise ValueError('Finite original FP32 components required')
        result=null_anchored_chain_set_loss(score,chains,labels)
        weighted=result['loss']*weights[ident]
        score['logits'].retain_grad(); weighted.backward()
        if any(float(score['logits'].grad[i])!=0. for i,y in enumerate(labels) if y==-1):
            raise ValueError('Unknown direct gain gradient must stay zero')
        modules=[('shared.'+k,getattr(student.shared,k)) for k in ('project','context','dna','fuse')]
        modules += [('head.'+k,getattr(student.head,k)) for k in A.HEAD_GROUPS]
        norms={}
        for name,module in modules:
            params=list(module.parameters())
            if name in ('head.junction','head.junction_additive') and not applicable:
                if any(p.grad is not None and not bool((p.grad==0).all()) for p in params):
                    raise ValueError('No known junction must give zero junction gradient')
                norms[name]=0.; continue
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual gradient: '+name)
            norms[name]=math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if norms[name]<=0: raise ValueError('Zero applicable branch gradient: '+name)
        exon_norm=float(student.head.exon[0].weight.grad[:,-3:].double().norm())
        junction_grad=student.head.junction[0].weight.grad
        junction_norm=float(junction_grad[:,-5:].double().norm()) if junction_grad is not None else 0.
        if exon_norm<=0 or (applicable and junction_norm<=0) or (not applicable and junction_norm!=0.):
            raise ValueError('Actual role columns must match known exon/junction coverage')
        if logits.requires_grad or any(p.requires_grad or p.grad is not None for p in source.parameters()):
            raise ValueError('Frozen local source must receive no gradient')
        torch.cuda.synchronize()
        receipts.append(dict(**witness,weighted_loss=float(weighted.detach()),gradient_norms=norms,
            known_junction_branch_applicable=applicable,role_exon_column_gradient_norm=exon_norm,
            role_junction_column_gradient_norm=junction_norm,unknown_direct_gain_gradient_zero=True,
            max_original_A31_local_score_error=float((reconstructed-expected).abs().max()) if len(chains) else 0.,
            frozen_local_source_gradients_absent=True,load_forward_backward_seconds=time.perf_counter()-started))
        print(receipts[-1],flush=True)
        del cached,features,dna,logits,reconstructed,expected,score,result,weighted
    if (len(receipts)!=8 or not all(any(r['species']==s and r['known_junction_branch_applicable'] for r in receipts)
            for s in R.C.R.ALLOW) or any(not torch.equal(v,states[name][k])
            for name,model in (('student',student),('source',source)) for k,v in model.state_dict().items())):
        raise ValueError('Unchanged fresh/frozen states and original eight-window coverage required')
    R.C.R.write_summary(OUT,dict(completed=True,experiment='M28-A40-ROLE-EVIDENCE-CONTEXT-PREFLIGHT',
        protocol=PROTOCOL,configuration=CONFIG,windows=1491,original_draws=1536,natural_candidates=16321,
        natural_counts=R.totals(records,records),
        replayed_original_draw_counts=R.totals(records,[d['window_id'] for d in draws]),
        eligible_reference_pairs=references,stronger_baseline_counts=baseline['selected_counts_per_species'],
        original_natural_identity_labels_and_stronger_baseline_reconciled=True,
        fresh_weights_unchanged=True,frozen_local_weights_unchanged=True,frozen_local_gradients_absent=True,
        original_A31_local_source_scores_reconciled=True,actual_role_column_gradients_verified=True,
        both_species_known_junction_gradient_coverage=True,student_supervised_checkpoint_read=False,
        frozen_supervised_A3_checkpoint_read=True,original_TRAIN_catalog_metadata_read=True,
        optimizer_steps_this_run=0,checkpoints_written=False,reference_injected_training_candidates=0,
        generator_decoding=False,generator_updated=False,new_candidate_search=False,new_backbone_forward=False,
        new_R4_labelling=False,DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        records=receipts,seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update role-conditioned token readiness, not trained performance, matched-capacity causality, '
            'novelty, raw-FASTA gene F1 or independent generalization. Frozen supervised A3 local evidence is an explicit input.'))


if __name__=='__main__':
    torch.set_num_threads(2); main()
