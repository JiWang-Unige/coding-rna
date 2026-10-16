#!/usr/bin/env python3
"""Reference-conditioned TRAIN mechanism probe, not natural gene prediction."""
import json
import math
import time
from collections import Counter
import torch
from torch.nn import functional as F
import run_full_train as R
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.query_pointer import STAGES

LIMITS = dict(beam_width=4, completed_limit_per_query=4, max_exons=128)


def coordinate_rank(coordinates, scores, target):
    """One-based descending-score rank; equal scores prefer smaller coordinate."""
    values = scores.detach().cpu().tolist()
    if len(coordinates) != len(values) or not all(math.isfinite(x) for x in values):
        raise ValueError('Finite coordinate-aligned scores required')
    position = coordinates.index(target)
    score = values[position]
    rank = 1 + sum(value > score or (value == score and x < target)
                   for x, value in zip(coordinates, values))
    return rank, max(values) - score


@torch.inference_mode()
def probe_query(window, sequence, chain, query_index):
    """Force the reference prefix and rank each NEXT decision with A5 semantics.

    Local top-4 passage is necessary, not sufficient: other prefix histories,
    previously completed chains and the score-bound stop can still remove it.
    No existence score, teacher assignment or reference-injected prediction.
    """
    sequence = sequence.upper()
    window.validate_chain(chain)
    if len(sequence) != window.length:
        raise ValueError('DNA and feature lengths differ')
    if len(chain) > LIMITS['max_exons']:
        return dict(query_index=query_index, grammar_supported=False,
                    unsupported_reason='exon_cap', events=[])
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    orf, ends = ORFIndex(sequence), _end_events(sequence)
    state, carry, previous_end = window.initial[query_index], '', -1
    events, first_exclusion = [], None
    boundary_sum, eos_sum = 0., 0.
    for depth, (a, b) in enumerate(chain):
        stage = 'start' if depth == 0 else 'next_start'
        lower = 0 if depth == 0 else previous_end + Budget().min_intron_bp
        entries = starts if depth == 0 else [x for x in acceptors if x >= lower]
        kind = 'stop' if depth == len(chain)-1 else 'donor'
        if a not in entries:
            return dict(query_index=query_index, grammar_supported=False,
                        unsupported_reason=stage, unsupported_exon=depth+1, events=events)
        entry_logp = F.log_softmax(window.pointer_logits(state, stage, lower), dim=-1)
        entry_rank, entry_deficit = coordinate_rank(entries, entry_logp[entries], a)
        geometric = list(range(lower, window.length))
        entry_geometric_rank, _ = coordinate_rank(geometric, entry_logp[geometric], a)
        state = window.advance(state, stage, a)
        choices = _legal_ends(orf, ends, a, carry, first_exon=depth == 0)
        target_indices = [k for k, (x, event_kind, _) in enumerate(choices)
                          if (x, event_kind) == (b, kind)]
        if len(target_indices) != 1:
            return dict(query_index=query_index, grammar_supported=False,
                        unsupported_reason='end_'+kind, unsupported_exon=depth+1, events=events)
        target = target_indices[0]
        coordinates = torch.tensor([x for x, _, _ in choices], device=state.device)
        tokens = (window.end_features[coordinates] + window.positions[coordinates] +
                  window.model.event.weight[STAGES['end']])
        states = window.model.recurrent(tokens, state.expand(len(choices), -1))
        eos = window.model.terminate(states).squeeze(-1)
        terminal = torch.tensor([k == 'stop' for _, k, _ in choices], device=state.device)
        branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
        end_logp = F.log_softmax(window.pointer_logits(state, 'end', a), dim=-1)
        scores = end_logp[coordinates] + branch
        same_kind = [k for k, (_, event_kind, _) in enumerate(choices) if event_kind == kind]
        end_rank, end_deficit = coordinate_rank(
            [choices[k][0] for k in same_kind], scores[same_kind], b)
        geometric = list(range(a+1, window.length+1))
        end_geometric_rank, _ = coordinate_rank(geometric, end_logp[geometric], b)
        entry_score, end_score, eos_score = float(entry_logp[a]), float(end_logp[b]), float(branch[target])
        boundary_sum += entry_score + end_score
        eos_sum += eos_score
        if first_exclusion is None:
            if entry_rank > LIMITS['beam_width']:
                first_exclusion = dict(exon=depth+1, stage=stage, rank=entry_rank)
            elif end_rank > (LIMITS['completed_limit_per_query'] if kind == 'stop' else LIMITS['beam_width']):
                first_exclusion = dict(exon=depth+1, stage='end_'+kind, rank=end_rank)
        events.append(dict(exon=depth+1, start=a, end=b, end_kind=kind,
                           entry_legal_rank=entry_rank, entry_geometric_rank=entry_geometric_rank,
                           entry_candidates=len(entries), entry_score_deficit=entry_deficit,
                           end_same_kind_rank=end_rank, end_geometric_rank=end_geometric_rank,
                           end_same_kind_candidates=len(same_kind), end_score_deficit=end_deficit,
                           entry_logp=entry_score, end_logp=end_score, eos_logp=eos_score,
                           prospective_eos_logit=float(eos[target]),
                           carry_before=carry, carry_after=choices[target][2]))
        state, carry, previous_end = states[target], choices[target][2], b
    return dict(query_index=query_index, grammar_supported=True, events=events,
                local_width4_necessary_conditions_passed=first_exclusion is None,
                first_local_exclusion=first_exclusion,
                boundary_log_score=boundary_sum, eos_log_score=eos_sum,
                path_log_score=boundary_sum+eos_sum, query_logit=float(window.query_logits[query_index]))


def aggregate(rows):
    return dict(reference_window_pairs=len(rows),
                natural_states=dict(Counter(r['state'] for r in rows)),
                grammar_supported=sum(r['grammar_supported'] for r in rows),
                any_query_passes_local_conditions=sum(bool(r['local_survivor_queries']) for r in rows),
                beam_missing_with_local_survivor=sum(not r['beam_free_hit'] and bool(r['local_survivor_queries']) for r in rows),
                beam_missing_all_queries_locally_excluded=sum(not r['beam_free_hit'] and r['grammar_supported'] and
                                                             not r['local_survivor_queries'] for r in rows))


@torch.inference_mode()
def main():
    import run_query_known_mean as A
    import analyze_query_beam_train as B
    source, out = B.OUT, R.ROOT/'outputs/M28-A5-TRUE-PREFIX-DIAGNOSTIC'
    original = json.loads((source/'summary.json').read_text())
    if ((source/'STATUS').read_text().strip() != 'COMPLETED' or not original['completed'] or
        original['limits'] != LIMITS or original['windows'] != 8 or
        original['pooled']['reference_window_pairs'] != 27 or original['DEV_forward_used'] or
        original['reference_used_during_decoding'] or original['test_setaria_used']):
        raise ValueError('Only the completed eight-window natural A5 TRAIN source is allowed')
    fitted = json.loads((A.STAGE_DIRS['fit']/'summary.json').read_text())
    A.check_last_identity(fitted)
    if original['checkpoint'] != fitted['final_checkpoint']:
        raise ValueError('A5 and fixed A3 last-weight identities differ')
    state = torch.load(R.ROOT/original['checkpoint'], map_location='cpu', weights_only=True)
    A.check_last_identity(state, state=True)
    model = A.new_model().eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True)
    del state
    predictions, references = R.rows(source/'predictions.jsonl'), R.rows(source/'references.jsonl')
    witnesses = R.rows(R.READY/'witnesses.jsonl')
    if (len(predictions) != 8 or len({r['window']['id'] for r in predictions}) != 8 or
        {r['window']['id'] for r in predictions} != {r['window_id'] for r in witnesses} or
        Counter(r['window']['species'] for r in predictions) != Counter({s:4 for s in R.ALLOW}) or
        len(references) != 27 or len({(r['window_id'], r['reference_id']) for r in references}) != 27):
        raise ValueError('Keep the same TRAIN windows and all 27 reference/window identities')
    selected = {r['window']['id'] for r in predictions}
    plan = {r['window']['id']:r for r in R.rows(R.PLAN/'feature_plan.jsonl') if r['window']['id'] in selected}
    began = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
    result, found_score_differences = [], []
    with (out/'query_profiles.jsonl').open('x') as log:
        for pred in predictions:
            w = pred['window']
            entry = plan[w['id']]
            refs = [r for r in references if r['window_id'] == w['id']]
            if (w['split'] != 'train' or w['seqid'] not in R.ALLOW[w['species']] or w['valid_bases'] != 24576 or
                any(w[k] != entry['window'][k] for k in R.GEOMETRY) or
                len(refs) != pred['expected_positive_chains'] or
                sum(r['exons'] for r in refs) != pred['expected_positive_exons'] or
                pred['reference_used_in_decoding']):
                raise ValueError('Unchanged TRAIN geometry and reference denominator required')
            cached = torch.load(R.ROOT/entry['artifact'], map_location='cpu', weights_only=True)
            if cached['sequence'] != pred['sequence']:
                raise ValueError('Source predictions and current TRAIN DNA differ')
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]
            window = model.query.for_window(hidden)
            greedy, beam = pred['greedy_queries'], pred['beam_queries']
            if [r['query_index'] for r in greedy] != list(range(16)) or [r['query_index'] for r in beam] != list(range(16)):
                raise ValueError('Frozen Q16 natural source changed')
            for ref in refs:
                chain = B.chain_key(ref['chain'])
                gh = [q['query_index'] for q in greedy if q['terminated'] and B.chain_key(q['chain']) == chain]
                bh = [q['query_index'] for q in beam if any(B.chain_key(c['chain']) == chain for c in q['completed'])]
                if (bool(gh), bool(bh)) != (ref['greedy_free_hit'], ref['beam_free_hit']):
                    raise ValueError('Natural hit identities changed')
                profiles = [probe_query(window, pred['sequence'], chain, q) for q in range(16)]
                for profile in profiles:
                    log.write(json.dumps(dict(window_id=w['id'], reference_id=ref['reference_id'],
                                              natural_greedy_hit=profile['query_index'] in gh,
                                              natural_beam_hit=profile['query_index'] in bh,
                                              reference_conditioned=True, **profile), allow_nan=False)+'\n')
                supported = all(p['grammar_supported'] for p in profiles)
                survivors = [p['query_index'] for p in profiles if p.get('local_width4_necessary_conditions_passed')]
                if not set(bh).issubset(survivors):
                    raise ValueError('A naturally returned query failed its necessary local ranks')
                for q in bh:
                    natural_score = next(c['path_log_score'] for c in beam[q]['completed'] if B.chain_key(c['chain']) == chain)
                    found_score_differences.append(abs(natural_score-profiles[q]['path_log_score']))
                result.append(dict(**ref, grammar_supported=supported, local_survivor_queries=survivors,
                                   natural_greedy_queries=gh, natural_beam_queries=bh,
                                   first_local_exclusion_by_query={str(p['query_index']):p.get('first_local_exclusion') for p in profiles},
                                   unsupported_reasons=sorted({p['unsupported_reason'] for p in profiles if not p['grammar_supported']})))
            log.flush()
            print(json.dumps(dict(window_id=w['id'], reference_pairs=len(refs), seconds=time.perf_counter()-began)), flush=True)
            del cached, features, dna, hidden, window
    if max(found_score_differences, default=0.) > .001:
        raise ValueError('Forced returned-path scores differ from the natural A5 source by more than 1e-3')
    with (out/'references.jsonl').open('x') as log:
        for row in result:
            log.write(json.dumps(row, allow_nan=False)+'\n')
    report = dict(experiment='M28-A5-TRUE-PREFIX-DIAGNOSTIC', completed=True,
                  source=str(source.relative_to(R.ROOT)), checkpoint=original['checkpoint'], limits=LIMITS,
                  windows=8, query_profiles=len(result)*16, pooled=aggregate(result),
                  per_species={s:aggregate([r for r in result if r['species'] == s]) for s in R.ALLOW},
                  by_exon_class={c:aggregate([r for r in result if r['exon_class'] == c]) for c in ('one','two_to_four','five_or_more')},
                  by_natural_state={s:aggregate([r for r in result if r['state'] == s]) for s in B.STATES},
                  returned_query_score_max_absolute_difference=max(found_score_differences, default=0.),
                  reference_conditioned=True, natural_predictions_reused_without_change=True,
                  optimizer_steps_this_run=0, checkpoints_written=False, new_backbone_forward=False,
                  DEV_forward_used=False, test_setaria_used=False,
                  seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
                  limitations='TRAIN pressure examples, true prefixes forced. Local passage is necessary not sufficient; not accuracy, causal attribution, a new decoder, F1/background or gate.')
    R.write_summary(out, report)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
