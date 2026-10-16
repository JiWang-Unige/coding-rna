#!/usr/bin/env python3
"""Read-only decomposition of the original A33 ordinal548 pointer replay failure."""
import json
import time
import torch
from torch.nn import functional as F
import run_emission_query_beam_train as A
from src.m28.query_pointer import STAGES

OUT = A.ROOT/'outputs/M28-A33-POINTER-NUMERICAL-DIAGNOSTIC'
ORDINAL = 548


def trace(window, sequence, row, batched):
    state = window.initial[row['query_index']]; carry = ''; previous = -1
    orf, ends = A.ORFIndex(sequence), A._end_events(sequence)
    grouped, components = 0., 0.; steps = []; states_saved = []
    for depth, (a, b) in enumerate(row['chain']):
        stage = 'start' if depth == 0 else 'next_start'
        lower = 0 if depth == 0 else previous+A.Budget().min_intron_bp
        lp = F.log_softmax(window.pointer_logits(state, stage, lower), -1)[a]
        state = window.advance(state, stage, a)
        choices = A._legal_ends(orf, ends, a, carry, depth == 0)
        kind = 'stop' if depth == len(row['chain'])-1 else 'donor'
        matches = [i for i, (x, k, _) in enumerate(choices) if (x, k) == (b, kind)]
        if len(matches) != 1:
            raise ValueError('Original observed branch outside raw grammar')
        k = matches[0]
        end_lp = F.log_softmax(window.pointer_logits(state, 'end', a), -1)[b]
        if batched:
            coordinates = torch.tensor([x for x, _, _ in choices], device=state.device)
            tokens = window.end_features[coordinates]+window.positions[coordinates]+window.model.event.weight[STAGES['end']]
            states = window.model.recurrent(tokens, state.expand(len(choices), -1))
            eos = window.model.terminate(states).squeeze(-1)
            terminal = torch.tensor([event == 'stop' for _, event, _ in choices], device=state.device)
            branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))[k]
            state = states[k]
        else:
            state = window.advance(state, 'end', b)
            eos = window.model.terminate(state).squeeze(-1)
            branch = F.logsigmoid(eos if kind == 'stop' else -eos)
        grouped += float(lp)+float(end_lp+branch)
        components += float(lp)+float(end_lp)+float(branch)
        steps.append(dict(exon=[a, b], event=kind, legal_ends=len(choices),
                          start_logp=float(lp), end_logp=float(end_lp), EOS_logp=float(branch)))
        states_saved.append(state.detach().clone()); carry, previous = choices[k][2], b
    return grouped, components, steps, states_saved


@torch.inference_mode()
def main():
    _, ids, index, pools, saved, _ = A.sources()
    w = pools[ORDINAL]['window']; ident = ids[ORDINAL]
    cached = torch.load(A.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
    features, dna = A.Q.window_inputs(cached, w, 'cuda'); model = A.C.load_generator()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    hidden = model(features, dna, [w['valid_bases']])['features'][0]; logits = model.segmentation(hidden)
    old = A.R.local_path_scores(logits, pools[ORDINAL]['chains'], w['valid_bases']).cpu()
    torch.testing.assert_close(old, torch.tensor(saved[ORDINAL]['variants']['A31_A3']['scores'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    queries = A.decode_emission_query_beam(model.query, hidden, logits, cached['sequence'], 128, 4, 4, w['valid_bases'])
    emitted = [r for q in queries for r in q['completed']]
    A.raw_grammar([r['chain'] for r in emitted], cached['sequence'])
    actual = A.R.local_path_scores(logits, [r['chain'] for r in emitted], w['valid_bases']).cpu()
    torch.testing.assert_close(actual, torch.tensor([r['emission_gain'] for r in emitted], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    window = model.query.for_window(hidden, w['valid_bases']); records = []
    for row in emitted:
        replay = A.replay_node_score(window, cached['sequence'], row['query_index'],
                                    dict(chain=row['chain'], entry=None, event=A.STOP))
        bg, bc, bs, batch_states = trace(window, cached['sequence'], row, True)
        sg, sc, ss, scalar_states = trace(window, cached['sequence'], row, False)
        records.append(dict(query_index=row['query_index'], rank=row['candidate_rank'], chain=row['chain'],
            native_pointer_score=row['path_log_score'], original_scalar_FP32_replay=float(replay),
            original_tolerance_pass=bool(torch.isclose(replay, replay.new_tensor(row['path_log_score']), rtol=1e-6, atol=1e-5)),
            batch_grouped_FP64=bg, scalar_grouped_FP64=sg,
            batch_components_FP64=bc, scalar_components_FP64=sc,
            batch_vs_native_error=bg-row['path_log_score'],
            scalar_vs_batch_grouped_error=sg-bg, scalar_FP32_vs_components_FP64_error=float(replay)-sc,
            max_end_state_batch_scalar_difference=max(float((x-y).abs().max()) for x, y in zip(batch_states, scalar_states)),
            batch_steps=bs, scalar_steps=ss))
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen generator changed or gradients appeared')
    with (OUT/'paths.jsonl').open('x') as handle:
        for row in records:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    summary = dict(completed=True, experiment='M28-A33-POINTER-NUMERICAL-DIAGNOSTIC',
        original_failed_job='13403217_0', ordinal=ORDINAL, window_id=ident, paths=len(records),
        original_tolerance_failures=sum(not r['original_tolerance_pass'] for r in records),
        max_batch_vs_native_error=max(abs(r['batch_vs_native_error']) for r in records),
        max_scalar_vs_batch_grouped_error=max(abs(r['scalar_vs_batch_grouped_error']) for r in records),
        max_scalar_FP32_accumulation_error=max(abs(r['scalar_FP32_vs_components_FP64_error']) for r in records),
        max_end_state_batch_scalar_difference=max(r['max_end_state_batch_scalar_difference'] for r in records),
        original_A31_scores_reproduced=True, raw_grammar_verified=True,
        final_emission_gains_directly_reproduced=True, model_weights_unchanged=True,
        optimizer_steps=0, new_parameters=0, new_backbone_forward=False, reference_or_catalog_opened=False,
        DEV_used=False, test_setaria_used=False, original_search_validator_or_protocol_changed=False,
        original_full_readout_restarted=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='One failed TRAIN window, numerical diagnosis only; no new effect labels or gate and no original run rescue.')
    A.C.R.write_summary(OUT, summary)


if __name__ == '__main__':
    torch.set_num_threads(2); main()
