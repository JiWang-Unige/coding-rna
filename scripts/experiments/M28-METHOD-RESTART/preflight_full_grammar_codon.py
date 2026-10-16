#!/usr/bin/env python3
"""Two fixed fresh TRAIN views, independent disposable source4608 updates."""
import argparse
import copy
import json
import os
import time
import torch
import fit_full_grammar_fresh_joint24 as B
import annotate_full_grammar_codon as A
from src.m28.full_grammar_codon import codon_coarse_channel_nll, centered_codon_scores
from src.m28.full_grammar_codon_skyline import decode_full_grammar_codon

C, J, G, Q, ROOT = B.C, B.J, B.G, B.Q, B.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-PREFLIGHT'
CONFIG = dict(B.CONFIG, **A.IDENTITY)


def checkpoint_fixture(state):
    A.validate_checkpoint(state)
    for key, value in (('codon_order', list(reversed(A.CODONS))), ('codon_centering', 'none'),
                       ('sampling_seed', 2), ('encoder_blocks_bp', 6144)):
        changed = copy.deepcopy(state)
        changed['configuration'][key] = value
        try:
            A.validate_checkpoint(changed)
        except ValueError:
            pass
        else:
            raise ValueError('Raw codon adapter accepted changed ' + key)
    changed = dict(state, raw_codon_scores=torch.zeros(63))
    try:
        A.validate_checkpoint(changed)
    except ValueError:
        pass
    else:
        raise ValueError('Raw codon adapter accepted wrong parameter shape')


def main(checks_job):
    began = time.perf_counter()
    checks = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-CHECKS'/checks_job
    if (checks/'STATUS').read_text().strip() != 'COMPLETED' or (checks/'JOBID').read_text().strip() != checks_job:
        raise ValueError('Completed independent codon enumeration checks required')
    plan, windows, catalogs = B.plan_and_population()
    selected = [next(i for i, r in enumerate(plan[:1536]) if r['species'] == species) for species in Q.ALLOW]
    C.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=A.PROTOCOL, checks_job=checks_job, selected_additional_steps=[i+1 for i in selected],
        selected_rule='same_first_original_fresh_draw_per_species_as_B_preflight'))
    source = C.starting_state()
    sequences, descriptors = C.P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    records = []
    for ordinal in selected:
        started = time.perf_counter()
        expected, w = plan[ordinal], windows[plan[ordinal]['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
        model, optimizer = C.load_model_optimizer(source)
        raw_scores = torch.nn.Parameter(torch.zeros(64, dtype=torch.float32, device='cpu'))
        optimizer.add_param_group({'params': [raw_scores]})
        joined = J.potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
        cuts = tuple(x for x in (offset, offset+w['valid_bases']) if 0 < x < len(context))
        result = codon_coarse_channel_nll(context, joined, 128, raw_scores, offset,
            offset+w['valid_bases'], targets, .1, cuts=cuts, checkpoint_bp=256)
        raw = result['loss']
        parent, _ = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
        if float(raw.detach()) != float(parent.detach()):
            raise ValueError('Zero codon scores changed actual source4608 native NLL')
        gradient, coefficient_gradient = torch.autograd.grad(raw, (joined, raw_scores), retain_graph=True)
        parent_gradient, = torch.autograd.grad(parent, joined, retain_graph=True)
        torch.testing.assert_close(gradient, parent_gradient, rtol=1e-6, atol=1e-7)
        if not bool(torch.isfinite(coefficient_gradient).all()):
            raise ValueError('Nonfinite actual codon VJP')
        total = raw * expected['objective_multiplier']
        if not bool(torch.isfinite(total)) or float(raw.detach()) < -1e-8:
            raise ValueError('Nonfinite or negative codon likelihood')
        total.backward()
        parameters = list(model.parameters()) + [raw_scores]
        if any(p.grad is None for p in parameters):
            raise ValueError('Missing native head or codon parameter gradient')
        norm = torch.nn.utils.clip_grad_norm_(parameters, 1., error_if_nonfinite=True)
        if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
            raise ValueError('Native pretrained encoder must remain frozen')
        optimizer.step()
        torch.cuda.synchronize()
        if (any(int(optimizer.state[p]['step']) != 4609 for p in model.parameters())
                or int(optimizer.state[raw_scores]['step']) != 1
                or any(not bool(torch.isfinite(p).all()) for p in parameters)
                or all(torch.equal(v.cpu(), source['model'][k]) for k, v in model.state_dict().items())):
            raise ValueError('Original head4609 and new codon1 finite disposable update required')
        fixture = dict(configuration=CONFIG, protocol=A.PROTOCOL, step=9216, epoch=6,
            encoder_revision=Q.F.REVISION, model=model.state_dict(), raw_codon_scores=raw_scores.detach().clone())
        checkpoint_fixture(fixture)
        restored_model, restored_scores = A.model_from_state(fixture, 'cuda')
        model.eval()
        with torch.inference_mode():
            trained = J.potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
            restored = J.raw_view(restored_model, backbone, tokenizer, chromosome, geometry, w['strand'])
        if not torch.equal(trained, restored) or not torch.equal(raw_scores.detach(), restored_scores):
            raise ValueError('Actual updated training/raw head or codon coefficients differ')
        left = decode_full_grammar_codon(context, trained, 128, raw_scores)
        right = decode_full_grammar_codon(context, restored, 128, restored_scores)
        for key in ('chains', 'gain', 'chain_gains'):
            if left[key] != right[key]:
                raise ValueError('Actual updated native training/raw codon MAP differs')
        txs, scopes = A.annotate_contig(restored_model, restored_scores, backbone, tokenizer,
            'toy14', 'ATGCGTNAGTTTAA', ordinal, 'cuda')
        if len(scopes) != 2 or any(s['windows'] != [(0, 14)] for s in scopes):
            raise ValueError('Actual native codon raw caller must retain both strands and all short DNA rows')
        row = dict(additional_step=ordinal+1, window_id=w['id'], species=w['species'], strand=w['strand'],
            native_context_bases=len(context), raw_NLL=float(raw.detach()), weighted_total=float(total.detach()),
            original_known_observation_denominator=expected['known_observation_denominator'],
            gradient_norm_before_clip=float(norm), initial_head_optimizer_step=4608,
            final_disposable_head_optimizer_step=4609, initial_codon_optimizer_step=0, final_codon_optimizer_step=1,
            zero_codon_native_NLL_parent_exact=True,
            zero_codon_role_VJP_parent_max_abs=float((gradient-parent_gradient).abs().max()),
            unnormalized_codon_VJP_l2=float(coefficient_gradient.double().norm()),
            updated_raw_codon_l2=float(raw_scores.detach().double().norm()),
            updated_centered_codon_l2=float(centered_codon_scores(raw_scores.detach()).norm()),
            actual_native_logits_raw_adapter_exact=True, actual_codon_scores_raw_adapter_exact=True,
            actual_updated_context_MAP_raw_adapter_exact=True, short_native_raw_both_strands_complete=True,
            short_complete_CDS_predictions=len(txs), seconds=time.perf_counter()-started)
        records.append(row)
        print(json.dumps(row, allow_nan=False), flush=True)
        del model, optimizer, joined, raw_scores, result, raw, total, gradient, parent_gradient, parent
        del restored_model, restored_scores, fixture, trained, restored, left, right, txs, scopes
    C.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, protocol=A.PROTOCOL, checks_job=checks_job, records=records,
        source_fit_job='13434603', support_job='13453408', original_TRAIN_DNA_sources=descriptors,
        disposable_updates_only=True, optimizer_steps_this_run=2, formal_training_started=False,
        checkpoints_written=False, raw_adapter_fixture_metadata_not_formal9216=True,
        saved_source_weights_and_optimizer_unchanged=True, pretrained_encoder_frozen_no_grad=True,
        cached_features_used=False, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began, interpretation='native_codon_energy_wiring_not_chain_accuracy'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.checks_job)
