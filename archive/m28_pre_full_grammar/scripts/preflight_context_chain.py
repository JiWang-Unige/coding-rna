#!/usr/bin/env python3
"""A37 natural-domain complete-chain supervision into fresh context, zero updates."""
from collections import Counter
import math
import time
import torch
import run_natural_chain_set as R
import readout_structured_emission_train as B
from diagnose_missing_path_emission_train import require_jobs
from src.m28.context_chain import ContextChain
from src.m28.chain_set import chain_set_loss

ROOT = R.ROOT
OUT = ROOT/'outputs/M28-A37-CONTEXT-CHAIN-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a37_context_chain_protocol.md'
CONFIG = dict(feature_dim=2048, hidden=128, seed=0, precision='FP32',
    model='fresh_SharedContext_jointly_learned_with_existing_JunctionChainHead',
    trainable_parameters=612324, context_parameters=364256, head_parameters=248068,
    initialization='one_seed0_stream_SharedContext_then_JunctionChainHead_no_old_weights',
    pretrained_encoder_frozen=True, source_natural_generator=R.C.GENERATOR,
    candidate_pool='original_16321_natural_A3_chains_no_reference_injection',
    loss='exact_compatible_known_candidate_set_NLL_divided_by_known_count',
    loss_unknown='exclude_-1_from_partition_and_gold_keep_in_natural_competition',
    score='exon_additive_plus_junction_additive_plus_chain_nonadditive_minus_window_null',
    threshold=0., threshold_scan=False, local_cross_entropy=False,
    draw_order_and_importance_weights='original_1536_draws_1491_unique_TRAIN',
    proposed_passes=3, proposed_optimizer_steps=4608,
    optimizer='AdamW', lr=3e-4, weight_decay=.01, betas=[.9,.999], eps=1e-8, clip_norm=1.,
    proposed_main_checkpoint='last_step_004608_only',
    comparison='frozen_A31_A3_same_natural_pool_and_same_zero_gain_operating_point')
HEAD_GROUPS = ('exon', 'sequence', 'additive', 'nonadditive', 'null', 'junction', 'junction_additive')


def new_model(device):
    torch.manual_seed(CONFIG['seed'])
    return ContextChain(CONFIG['feature_dim'], CONFIG['hidden']).to(device)


def population():
    draws, index, windows, records, source = R.population()
    labels = Counter(r['label'] for rec in records.values() for r in rec['candidates'])
    if (len(draws) != 1536 or len(records) != 1491 or labels != {1:2380, 0:12566, -1:1375} or
            any(r['origin'] != 'free' for rec in records.values() for r in rec['candidates'])):
        raise ValueError('Complete unchanged original natural-only TRAIN cohort required')
    require_jobs(['13398521', '13398524'])
    old = B.R.completed(B.R.READOUT, '13398521'); checked = B.R.completed(B.R.CHECK, '13398524')
    if (old['selected_counts_per_species'] != checked['selected_counts_per_species'] or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] or
            not checked['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed']):
        raise ValueError('Actual closed stronger A31_A3 baseline required')
    previous = R.C.R.rows(B.R.READOUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(B.R.READOUT/'windows.jsonl')
    if [row['window']['id'] for row in previous] != list(records) or len(labelled) != 1491:
        raise ValueError('Original stronger-baseline ordered cohort differs')
    references = Counter()
    for prior, lab, (ident, rec) in zip(previous, labelled, records.items()):
        if (lab['window'] != prior['window'] or prior['reference_used_in_scoring'] or
                any(rec['window'][k] != prior['window'][k] for k in R.C.R.GEOMETRY) or
                [r['chain'] for r in rec['candidates']] != prior['variants']['A31_A3']['chains'] or
                [r['label'] for r in rec['candidates']] != lab['labels']['A31_A3']):
            raise ValueError('Original natural identity and stronger-baseline label pairing differs')
        references.update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
    if references != {'eligible':2916, 'five_or_more':1086}:
        raise ValueError('Original full reference/window denominators differ')
    return draws, index, windows, records, source, old, dict(references)


def main():
    require_jobs(['13406144', '13408449', '13408451', '13409553'])
    closed = B.R.completed(B.CHECK, '13408451')
    if closed['TRAIN_preparation_supported'] is not False:
        raise ValueError('A36 closed FAIL must remain separate from this new model')
    draws, index, windows, records, source, baseline, references = population()
    model = new_model('cuda').train()
    if (sum(p.numel() for p in model.parameters()) != CONFIG['trainable_parameters'] or
            sum(p.numel() for p in model.shared.parameters()) != CONFIG['context_parameters'] or
            sum(p.numel() for p in model.head.parameters()) != CONFIG['head_parameters']):
        raise ValueError('Declared fresh joint-context model structure differs')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    weights = {d['window_id']: d['importance_weight_uniform_within_stratum'] for d in draws}
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats(); receipts = []
    for witness in R.P.select_gradient_windows(records):
        ident = witness['window_id']; w = windows[ident]; local = records[ident]['candidates']
        chains, labels = [r['chain'] for r in local], [r['label'] for r in local]
        applicable = any(y >= 0 and len(c) > 1 for c, y in zip(chains, labels))
        started = time.perf_counter()
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        features, dna = R.C.R.Q.window_inputs(cached, w, 'cuda')
        if features.requires_grad or dna.requires_grad:
            raise ValueError('Original cached encoder features and DNA must remain detached')
        model.zero_grad(set_to_none=True)
        score = model(features, dna, [w['valid_bases']], chains)
        if any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in score.values()):
            raise ValueError('Finite original FP32 joint-context scores required')
        result = chain_set_loss(score, chains, labels); weighted = result['loss']*weights[ident]
        if not bool(torch.isfinite(weighted)):
            raise ValueError('Nonfinite actual natural-domain structured loss')
        score['logits'].retain_grad(); weighted.backward()
        if any(float(score['logits'].grad[i]) != 0. for i, y in enumerate(labels) if y == -1):
            raise ValueError('Original unknown direct gain gradient must remain zero')
        modules = [('shared.'+k, getattr(model.shared, k)) for k in ('project','context','dna','fuse')]
        modules += [('head.'+k, getattr(model.head, k)) for k in HEAD_GROUPS]
        norms = {}
        for name, module in modules:
            params = list(module.parameters())
            if name in ('head.junction','head.junction_additive') and not applicable:
                if any(p.grad is not None and not bool((p.grad == 0).all()) for p in params):
                    raise ValueError('No known junction must give no junction update: '+name)
                norms[name] = 0.; continue
            if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                raise ValueError('Missing/nonfinite actual gradient: '+name)
            norms[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
            if norms[name] <= 0:
                raise ValueError('Zero actual applicable branch gradient: '+name)
        torch.cuda.synchronize()
        receipts.append(dict(**witness, weighted_loss=float(weighted.detach()), gradient_norms=norms,
            known_junction_branch_applicable=applicable, unknown_direct_gain_gradient_zero=True,
            load_forward_backward_seconds=time.perf_counter()-started,
            **{k: v for k, v in result.items() if not torch.is_tensor(v)}))
        print(receipts[-1], flush=True)
        del cached, features, dna, score, result, weighted
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            not all(any(r['species'] == s and r['known_junction_branch_applicable'] for r in receipts)
                    for s in R.C.R.ALLOW)):
        raise ValueError('Unchanged state and actual known junction coverage in both species required')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A37-CONTEXT-CHAIN-PREFLIGHT',
        protocol=PROTOCOL, configuration=CONFIG, windows=1491, original_draws=1536,
        natural_candidates=16321, natural_counts=R.totals(records, records),
        replayed_original_draw_counts=R.totals(records, [d['window_id'] for d in draws]),
        windows_losing_all_injected_gold=sum(bool(source[i]['counts']['label_1']) and
                                            not records[i]['counts']['label_1'] for i in records),
        eligible_reference_pairs=references, stronger_baseline_counts=baseline['selected_counts_per_species'],
        original_natural_identity_labels_and_stronger_baseline_reconciled=True,
        all_natural_gold_sets_compatible=True, fresh_weights_unchanged=True,
        both_species_known_junction_gradient_coverage=True, original_TRAIN_catalog_metadata_read=True,
        supervised_checkpoint_read=False, generator_forward=False, local_segmentation_head=False,
        optimizer_steps_this_run=0, checkpoints_written=False, new_candidate_search=False,
        new_backbone_forward=False, new_R4_labelling=False, DEV_used=False,
        test_setaria_used=False, A33_outputs_opened=False, reference_injected_training_candidates=0,
        records=receipts, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Zero-update joint-context gradient/identity readiness only; no trained performance, '
            'novelty, matched-capacity causal ablation or independent generalization.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
