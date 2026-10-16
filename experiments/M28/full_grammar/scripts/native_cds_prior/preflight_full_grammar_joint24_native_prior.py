#!/usr/bin/env python3
"""Original-TRAIN disposable updates for the joint24576 dual-native-prior composition."""
import argparse
import io
import json
import math
import os
import time
import torch
from torch import nn
import fit_full_grammar_joint24 as J
import fit_full_grammar_native_cds_prior as N
import annotate_full_grammar_joint24_native_cds_prior as A
from native_cds_role_prior import NativeCDSRolePrior

C = J.C
ROOT = J.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-PREFLIGHT'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-CHECKS'
CONFIG = dict(N.CONFIG,**A.IDENTITY,matched_controls=['STRUCTURED_CONTINUE_final9216',
    'NATIVE_CDS_PRIOR_6144_final9216','JOINT24576_C0_final9216'])


class ZeroPriorReplay(nn.Module):
    """Preflight only: compare the zero prior to C0 without a second native forward."""
    def __init__(self,adapter):
        super().__init__(); self.adapter = adapter; self.windows_replayed = 0

    def forward(self,features,dna,valid_bases):
        result = self.adapter(features,dna,valid_bases)
        with torch.no_grad():
            plain = C.G.M28Core.forward(self.adapter.core,features,dna,valid_bases)
        if not torch.equal(result['segmentation_logits'].detach(),plain['segmentation_logits']):
            raise ValueError('Zero-initialized combined prior must exactly preserve joint24576 C0 logits')
        self.windows_replayed += 1
        return result


def run(checks_job):
    began = time.perf_counter()
    checked = json.loads((CHECKS/'summary.json').read_text())
    source = json.loads((J.PREFLIGHT/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['job_id'] != checks_job or checked['identity'] != A.IDENTITY):
        raise ValueError('Completed matching combined-caller artificial checks required')
    if ((J.PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
            or source['job_id'] != '13453758' or source['configuration'] != J.CONFIG
            or not source['all_eight_joint24576_potentials_exact']):
        raise ValueError('Completed original joint24576 native preflight required')
    state = C.starting_state()
    draws,windows,catalogs,plan = J.original_plan()
    selected = [next(i for i,d in enumerate(draws) if d['species'] == s) for s in C.Q.ALLOW]
    original = {r['draw_ordinal']:r for r in source['records']}
    C.write(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=A.PROTOCOL,checks_job=checks_job,selected_draw_ordinals=selected))
    sequences,descriptors = C.P.read_train_DNA()
    full,tokenizer = A.N.load_native()
    records = []
    for ordinal in selected:
        started = time.perf_counter()
        w = windows[draws[ordinal]['window_id']]
        chromosome,geometry,context,offset,targets = C.example(w,plan[ordinal],sequences,catalogs)
        core,optimizer = N.load_model_optimizer(state)
        adapter = A.N.NativeCDSHeadAdapter(core,full.score)
        zero_replay = ZeroPriorReplay(adapter)
        joined = J.potentials(zero_replay,full.model,tokenizer,chromosome,geometry,w['strand'])
        raw,components = C.objective('STRUCTURED_CONTINUE',context,joined,offset,w['valid_bases'],targets)
        if not math.isclose(float(raw.detach()),original[ordinal]['raw_NLL'],rel_tol=1e-6,abs_tol=1e-9):
            raise ValueError('Zero prior starting loss differs from original joint24576 C0 preflight')
        total,norm = C.update(core,optimizer,full,raw,plan[ordinal]['objective_multiplier'])
        prior = core.native_cds_role_prior.weight
        old = [p for name,p in core.named_parameters() if name != 'native_cds_role_prior.weight']
        if (any(int(optimizer.state[p]['step']) != 4609 for p in old)
                or int(optimizer.state[prior]['step']) != 1 or not bool(prior.detach().count_nonzero())):
            raise ValueError('Old4608->4609 and new0->1 independent disposable update required')
        buffer = io.BytesIO(); torch.save(core.state_dict(),buffer); buffer.seek(0)
        restored = NativeCDSRolePrior().to('cuda').eval()
        restored.load_state_dict(torch.load(buffer,map_location='cuda',weights_only=True),strict=True)
        core.eval()
        restored_adapter = A.N.NativeCDSHeadAdapter(restored,full.score).eval()
        with torch.inference_mode():
            actual = J.potentials(adapter,full.model,tokenizer,chromosome,geometry,w['strand'])
            replay = J.raw_view(restored_adapter,full.model,tokenizer,chromosome,geometry,w['strand'])
        if not torch.equal(actual,replay):
            raise ValueError('Updated neural weights and raw stitching must exactly replay the same joint-native potentials')
        row = dict(draw_ordinal=ordinal,window_id=w['id'],species=w['species'],raw_NLL=float(raw.detach()),
            weighted_total=total,gradient_norm_before_clip=norm,old_optimizer_step=4609,
            new_prior_optimizer_step=1,prior_weights=prior.detach().cpu().tolist(),
            zero_prior_original_joint_window_logits_exact=zero_replay.windows_replayed,
            zero_prior_original_joint_loss_matches=True,strict_neural_weight_in_memory_roundtrip=True,
            updated_raw_adapter_potentials_exact=True,seconds=time.perf_counter()-started)
        records.append(row); print(json.dumps(row,allow_nan=False),flush=True)
        del raw,components,joined,zero_replay,optimizer,buffer,actual,replay,targets
    # Native short raw call uses the last disposable neural weights, not a formal checkpoint.
    with torch.inference_mode():
        txs,scopes = A.J.annotate_contig(restored_adapter,full.model,tokenizer,'toy14','ATGCGTNAGTTTAA',0,'cuda')
    if len(scopes) != 2 or any(s['windows'] != [(0,14)] for s in scopes):
        raise ValueError('Actual short raw call must retain both strands and all fourteen bases')
    if any(p.requires_grad or p.grad is not None for p in full.parameters()):
        raise ValueError('Native encoder and both heads must remain frozen without gradients')
    C.write(OUT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol=A.PROTOCOL,checks_job=checks_job,source_joint_preflight_job='13453758',records=records,
        native_short_raw=dict(input_bases=14,strands=2,complete_CDS_predictions=len(txs),reference_used=False),
        original_TRAIN_DNA_sources=descriptors,optimizer_steps_this_run=2,disposable_updates_only=True,
        checkpoints_written=False,formal_fit_started=False,saved_source_weights_and_optimizer_unchanged=True,
        pretrained_encoder_and_heads_frozen_no_grad=True,cached_features_used=False,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,interpretation='combined_native_wiring_not_chain_accuracy_or_interaction_benefit'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-job',required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    run(args.checks_job)
