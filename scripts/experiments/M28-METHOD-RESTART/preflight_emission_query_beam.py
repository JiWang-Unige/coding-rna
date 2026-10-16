#!/usr/bin/env python3
"""A33 zero-update search preflight on the original eight TRAIN metadata cases."""
import json
import time
import torch
import run_local_path_caller as R
from query_search_score_replay import replay_node_score, STOP
from src.m28.emission_query_beam import decode_emission_query_beam, augment_pool
from src.m28.local_path_score import local_path_scores
from src.m28.phase_evidence_chain import candidate_roles
from src.m28.core import select_nonoverlapping

ROOT, C, Q = R.ROOT, R.C, R.Q
OUT = ROOT/'outputs/M28-A33-EMISSION-BEAM-PREFLIGHT-R1'
PROTOCOL = 'reports/M28-METHOD-RESTART/a33_emission_beam_protocol.md'
CONFIG = dict(generator_checkpoint=C.GENERATOR, generator_optimizer_steps=4608,
    beam_width=4, completed_limit=4, max_exons=128, local_pointer_topk_unchanged=True,
    global_retention='A31_role_emission_gain_then_original_pointer_log_score',
    entry_preview='START3_or_closed_prefix_plus_intron_and_acceptor_first_base',
    stop_bound='all_retained_live_prefix_gain_plus_suffix_max_positive_role_gain',
    pool='original_A31_A3_pool_union_emission_beam4', all_role_weights=1., null=0.,
    model_precision='FP32', score_accumulation='FP64', new_parameters=0, optimizer_steps=0)


@torch.inference_mode()
def main():
    prepared = R.completed(R.READOUT, 13398521); checked = R.completed(R.CHECK, 13398524)
    if (prepared['configuration'] != R.CONFIG or checked['configuration'] != R.CONFIG or
            prepared['windows'] != 1491 or checked['windows'] != 1491 or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified']):
        raise ValueError('Closed original A31 TRAIN source required')
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original complete TRAIN metadata required')
    cases = Q.select_timing_draws(draws, windows, catalogs)
    if len(cases) != 8:
        raise ValueError('Original eight metadata cases required')
    pool_source = R.completed(C.FIT); C.validate_head(pool_source)
    pools = {r['window']['id']: r for r in C.R.rows(C.FIT/'free_pool.jsonl')}
    saved = {r['window']['id']: r for r in C.R.rows(R.READOUT/'natural_scores.jsonl')}
    if (set(pools) != set(windows) or set(saved) != set(windows) or
            sum(len(r['chains']) for r in pools.values()) != 16321):
        raise ValueError('Original whole TRAIN natural pool and A31 score population required')
    model = C.load_generator()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649457:
        raise ValueError('Same frozen A3 parameter population required')
    records = []; began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for case in cases:
            ident = case['draw']['window_id']; w = windows[ident]; free = pools[ident]
            if (w['split'] != 'train' or w['seqid'] not in C.R.ALLOW[w['species']] or
                    free['reference_used_in_pool_generation'] or
                    any(w[k] != free['window'][k] for k in C.R.GEOMETRY)):
                raise ValueError('Original TRAIN-only free geometry required')
            started = time.perf_counter()
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            Q.validate_cached(cached, w)
            features, dna = Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]
            logits = model.segmentation(hidden)
            if logits.dtype != torch.float32:
                raise ValueError('Frozen FP32 model required')
            old_chains, old_sources = C.free_pool(model, hidden, cached['sequence'])
            old_chains = [[list(exon) for exon in chain] for chain in old_chains]
            if old_chains != free['chains'] or old_sources != free['decoder_sources']:
                raise ValueError('Actual old natural geometry/order/sources must reproduce')
            baseline = saved[ident]['variants']['A31_A3']
            if baseline['chains'] != old_chains:
                raise ValueError('Actual old A31 scored pool differs')
            old_scores = local_path_scores(logits, old_chains, w['valid_bases'])
            torch.testing.assert_close(old_scores.cpu(), torch.tensor(baseline['scores'], dtype=torch.float64),
                                       rtol=1e-12, atol=1e-8)
            t = time.perf_counter()
            queries = decode_emission_query_beam(model.query, hidden, logits, cached['sequence'], 128, 4, 4, w['valid_bases'])
            torch.cuda.synchronize(); search_seconds = time.perf_counter()-t
            emitted = [r for q in queries for r in q['completed']]
            direct = local_path_scores(logits, [r['chain'] for r in emitted], w['valid_bases'])
            actual = torch.tensor([r['emission_gain'] for r in emitted], dtype=torch.float64)
            torch.testing.assert_close(actual, direct.cpu(), rtol=1e-12, atol=1e-8)
            # Independent original raw-grammar and pointer/EOS score replay.
            window = model.query.for_window(hidden, w['valid_bases']); replay_errors = []
            for row in emitted:
                node = dict(chain=row['chain'], entry=None, event=STOP)
                pointer = replay_node_score(window, cached['sequence'], row['query_index'], node)
                torch.testing.assert_close(pointer, pointer.new_tensor(row['path_log_score']), rtol=1e-6, atol=1e-5)
                replay_errors.append(abs(float(pointer)-row['path_log_score']))
            chains, sources = augment_pool(old_chains, old_sources, queries)
            scores = local_path_scores(logits, chains, w['valid_bases'])
            logp = logits[:w['valid_bases']].double().log_softmax(-1)
            direct_logp = []
            for chain in chains:
                roles = candidate_roles(chain, w['valid_bases'], logits.device)
                direct_logp.append(float((logp[torch.arange(w['valid_bases'], device=logits.device), roles]-logp[:, 0]).sum()))
            torch.testing.assert_close(scores.cpu(), torch.tensor(direct_logp, dtype=torch.float64), rtol=1e-12, atol=1e-8)
            if not {tuple(map(tuple, c)) for c in old_chains} <= {tuple(map(tuple, c)) for c in chains}:
                raise ValueError('Original candidate geometry must remain in augmented pool')
            value = dict(chains=chains, scores=scores.cpu().tolist(),
                         selected=select_nonoverlapping(chains, scores.cpu().tolist()))
            handle.write(json.dumps(dict(window=w, original_A31_A3=baseline, A33=value,
                decoder_sources=sources, emission_queries=queries,
                reference_used_in_pool_generation=False, reference_used_in_scoring=False), allow_nan=False)+'\n')
            torch.cuda.synchronize()
            record = dict(window_id=ident, species=w['species'], category=case['category'], draw_ordinal=case['ordinal'],
                valid_bases=w['valid_bases'], old_candidates=len(old_chains),
                emission_completed_query_paths=len(emitted), augmented_candidates=len(chains),
                newly_added_unique_chains=len(chains)-len(old_chains),
                search_seconds=search_seconds, full_case_seconds=time.perf_counter()-started,
                max_incremental_vs_direct_error=max((abs(float(x)-float(y)) for x, y in zip(actual, direct.cpu())), default=0.),
                max_direct_log_probability_error=max((abs(float(x)-y) for x, y in zip(scores.cpu(), direct_logp)), default=0.),
                max_pointer_replay_error=max(replay_errors, default=0.),
                query_counts=[q['counts'] for q in queries], original_natural_geometry_order_sources_reproduced=True,
                original_A31_local_scores_reproduced=True, returned_raw_grammar_verified=True)
            records.append(record)
            print(json.dumps({k: v for k, v in record.items() if k != 'query_counts'}, allow_nan=False), flush=True)
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen A3 weights changed or gradients appeared')
    C.R.write_summary(OUT, dict(completed=True, experiment='M28-A33-EMISSION-BEAM-PREFLIGHT',
        configuration=CONFIG, protocol=PROTOCOL, cases=8, records=records,
        original_TRAIN_draws=1536, unique_TRAIN_windows=1491, model_parameters=649457,
        original_natural_pool_and_A31_scores_reproduced=True, incremental_final_role_scores_verified=True,
        returned_raw_grammar_verified=True, original_pool_geometry_preserved=True,
        model_weights_unchanged=True, all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0, new_backbone_forward=False, feature_cache_used=True,
        new_reference_labelling=False, reference_injected_candidates=0, DEV_used=False, test_setaria_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began,
        limitation='Eight metadata cases, engineering and timing only. New chain count is not correct-chain '
            'coverage or an effect gate. Global retention changes, local pointer top-k does not. No training, '
            'full TRAIN gate, gene F1, native FASTA run, independent generalization or publication PASS.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
