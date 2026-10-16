#!/usr/bin/env python3
"""A16 zero-update legal-support mass diagnostic on the closed A5 TRAIN cohort."""
from collections import Counter
import json
import math
import time
import torch
from torch.nn import functional as F
import run_chain_quality as C
import run_full_train as R
import analyze_query_beam_train as B
import analyze_true_prefix_train as T
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.query_pointer import STAGES
from src.m28.legal_support import conditional_support

OUT = R.ROOT/'outputs/M28-A16-LEGAL-SUPPORT-DIAGNOSTIC'
SOURCE = R.ROOT/'outputs/M28-A5-TRUE-PREFIX-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a16_legal_support_protocol.md'


def support(log_weights, target):
    result = conditional_support(log_weights, target)
    values = {k: float(v) for k, v in result.items()}
    if not all(math.isfinite(v) for v in values.values()) or values['log_mass'] > 1e-5:
        raise ValueError('The actual legal event support must have finite mass at most one')
    return values


def probe(window, sequence, chain, query_index, original):
    window.validate_chain(chain); sequence = sequence.upper()
    starts = [i for i in range(len(sequence)-2) if sequence[i:i+3] == 'ATG']
    acceptors = [i+2 for i in range(len(sequence)-2) if sequence[i:i+2] == 'AG']
    ends, orf = _end_events(sequence), ORFIndex(sequence)
    state, carry, previous = window.initial[query_index], '', -1
    events = []; raw_total = 0.; conditional_total = 0.; maximum_delta = 0.
    for depth, (a, b) in enumerate(chain):
        stage = 'start' if depth == 0 else 'next_start'
        lower = 0 if depth == 0 else previous+Budget().min_intron_bp
        teacher_lower = 0 if depth == 0 else previous+1
        entries = starts if depth == 0 else [x for x in acceptors if x >= lower]
        if a not in entries:
            raise ValueError('Original A5 target lost decoder grammar support')
        entry_logp = F.log_softmax(window.pointer_logits(state, stage, lower), dim=-1)
        teacher_logp = F.log_softmax(window.pointer_logits(state, stage, teacher_lower), dim=-1)
        entry = support(entry_logp[entries], entries.index(a))
        entry_rank, _ = T.coordinate_rank(entries, entry_logp[entries], a)
        state = window.advance(state, stage, a)
        choices = _legal_ends(orf, ends, a, carry, first_exon=depth == 0)
        kind = 'stop' if depth == len(chain)-1 else 'donor'
        targets = [i for i, (x, k, _) in enumerate(choices) if (x, k) == (b, kind)]
        if len(targets) != 1:
            raise ValueError('Original A5 end/carry support changed')
        target = targets[0]
        coordinates = torch.tensor([x for x, _, _ in choices], device=state.device)
        tokens = window.end_features[coordinates]+window.positions[coordinates]+window.model.event.weight[STAGES['end']]
        states = window.model.recurrent(tokens, state.expand(len(choices), -1))
        eos = window.model.terminate(states).squeeze(-1)
        terminal = torch.tensor([k == 'stop' for _, k, _ in choices], device=state.device)
        branch = torch.where(terminal, F.logsigmoid(eos), F.logsigmoid(-eos))
        end_logp = F.log_softmax(window.pointer_logits(state, 'end', a), dim=-1)
        joint = end_logp[coordinates]+branch
        end = support(joint, target)
        end_boundary_log_mass = float(torch.logsumexp(end_logp[coordinates.unique()].double(), 0))
        same_kind = [i for i, (_, k, _) in enumerate(choices) if k == kind]
        end_rank, _ = T.coordinate_rank([choices[i][0] for i in same_kind], joint[same_kind], b)
        old = original['events'][depth]
        raw = dict(entry_logp=float(entry_logp[a]), end_logp=float(end_logp[b]), eos_logp=float(branch[target]))
        if (old['exon'], old['start'], old['end'], old['end_kind'], old['carry_before'], old['carry_after'],
            old['entry_legal_rank'], old['end_same_kind_rank']) != (
            depth+1, a, b, kind, carry, choices[target][2], entry_rank, end_rank):
            raise ValueError('A5 true-prefix geometry, ranks or real carry changed')
        maximum_delta = max(maximum_delta, *(abs(raw[k]-old[k]) for k in raw))
        raw_step = sum(raw.values())
        conditional_step = entry['conditional_target_logp']+end['conditional_target_logp']
        if abs(conditional_step-(raw_step-entry['log_mass']-end['log_mass'])) > 1e-5:
            raise ValueError('Raw/conditional score decomposition differs')
        raw_total += raw_step; conditional_total += conditional_step
        events.append(dict(exon=depth+1, stage=stage, start=a, end=b, end_kind=kind,
            entry_legal_candidates=len(entries), end_legal_events=len(choices),
            entry_geometry_positions=window.length-lower, end_geometry_positions=window.length-a,
            entry_log_mass=entry['log_mass'], end_boundary_log_mass=end_boundary_log_mass,
            end_joint_log_mass=end['log_mass'],
            entry_conditional_target_logp=entry['conditional_target_logp'],
            end_joint_conditional_target_logp=end['conditional_target_logp'],
            teacher_entry_target_logp=float(teacher_logp[a]), decoder_entry_target_logp=raw['entry_logp'],
            decoder_vs_teacher_lower_target_logp_change=raw['entry_logp']-float(teacher_logp[a]),
            entry_legal_rank=entry_rank, end_same_kind_rank=end_rank,
            carry_before=carry, carry_after=choices[target][2], **raw))
        state, carry, previous = states[target], choices[target][2], b
    maximum_delta = max(maximum_delta, abs(raw_total-original['path_log_score']))
    if maximum_delta > 1e-3:
        raise ValueError('A5 unchanged scores failed replay')
    return dict(query_index=query_index, reference_conditioned=True, events=events,
        raw_path_log_score=raw_total, hypothetical_conditioned_path_log_score=conditional_total,
        original_score_max_absolute_difference=maximum_delta,
        local_width4_necessary_conditions_passed=original['local_width4_necessary_conditions_passed'])


def reference_means(profiles):
    events = [e for p in profiles for e in p['events']]
    groups = {'start': [e for e in events if e['stage'] == 'start'],
              'next_start': [e for e in events if e['stage'] == 'next_start'],
              'end': events}
    return {stage: dict(events=len(values),
        mean_log_mass=sum(e['end_joint_log_mass' if stage == 'end' else 'entry_log_mass'] for e in values)/len(values),
        mean_probability_mass=sum(math.exp(e['end_joint_log_mass' if stage == 'end' else 'entry_log_mass']) for e in values)/len(values),
        mean_decoder_teacher_lower_shift=(0. if stage == 'end' else
            sum(e['decoder_vs_teacher_lower_target_logp_change'] for e in values)/len(values)))
        if values else None for stage, values in groups.items()}


def aggregate(rows):
    result = dict(reference_window_pairs=len(rows))
    for stage in ('start', 'next_start', 'end'):
        values = [r['support'][stage] for r in rows if r['support'][stage] is not None]
        result[stage] = dict(reference_pairs_with_events=len(values),
            **{k: sum(v[k] for v in values)/len(values) if values else None for k in
               ('mean_log_mass', 'mean_probability_mass', 'mean_decoder_teacher_lower_shift')})
    return result


@torch.inference_mode()
def main():
    source = json.loads((SOURCE/'summary.json').read_text())
    if ((SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or (SOURCE/'JOBID').read_text().strip() != '13359200' or
        not source['completed'] or source['checkpoint'] != C.GENERATOR or source['windows'] != 8 or
        source['query_profiles'] != 432 or source['pooled']['reference_window_pairs'] != 27 or
        source['DEV_forward_used'] or source['test_setaria_used']):
        raise ValueError('Only the closed A5 true-prefix TRAIN source is allowed')
    predictions = R.rows(B.OUT/'predictions.jsonl'); references = R.rows(SOURCE/'references.jsonl')
    old = R.rows(SOURCE/'query_profiles.jsonl')
    keyed = {(r['window_id'], r['reference_id'], r['query_index']): r for r in old}
    if (len(predictions) != 8 or len({r['window']['id'] for r in predictions}) != 8 or
        Counter(r['window']['species'] for r in predictions) != Counter({s: 4 for s in R.ALLOW}) or
        len(references) != 27 or len(old) != 432 or len(keyed) != 432):
        raise ValueError('Original eight windows/27 references/432 query profiles required')
    expected = {(r['window_id'], r['reference_id'], q) for r in references for q in range(16)}
    if set(keyed) != expected or any(not r['grammar_supported'] for r in old):
        raise ValueError('Keep every original supported reference/query, without performance filtering')
    selected = {r['window']['id'] for r in predictions}
    if selected != {r['window_id'] for r in R.rows(R.READY/'witnesses.jsonl')}:
        raise ValueError('The original metadata-selected eight-window set changed')
    plan = {r['window']['id']: r for r in R.rows(R.PLAN/'feature_plan.jsonl') if r['window']['id'] in selected}
    generator = C.load_generator(); began = time.perf_counter(); results = []; maximum_delta = 0.
    torch.cuda.reset_peak_memory_stats(); profile_count = 0; event_count = 0
    with (OUT/'query_profiles.jsonl').open('x') as log:
        for pred in predictions:
            w = pred['window']; entry = plan[w['id']]
            if (w['split'] != 'train' or w['seqid'] not in R.ALLOW[w['species']] or pred['reference_used_in_decoding'] or
                any(w[k] != entry['window'][k] for k in R.GEOMETRY)):
                raise ValueError('Only unchanged original TRAIN geometry')
            cached = torch.load(R.ROOT/entry['artifact'], map_location='cpu', weights_only=True)
            if cached['sequence'] != pred['sequence']:
                raise ValueError('Original oriented TRAIN DNA changed')
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            hidden = generator(features, dna, [w['valid_bases']])['features'][0]
            window = generator.query.for_window(hidden, w['valid_bases'])
            local_refs = [r for r in references if r['window_id'] == w['id']]
            if len(local_refs) != pred['expected_positive_chains']:
                raise ValueError('Original reference denominator changed')
            for ref in local_refs:
                profiles = [probe(window, pred['sequence'], ref['chain'], q,
                                  keyed[(w['id'], ref['reference_id'], q)]) for q in range(16)]
                for profile in profiles:
                    log.write(json.dumps(dict(window_id=w['id'], reference_id=ref['reference_id'], **profile), allow_nan=False)+'\n')
                    maximum_delta = max(maximum_delta, profile['original_score_max_absolute_difference'])
                    profile_count += 1; event_count += len(profile['events'])
                results.append(dict(**ref, support=reference_means(profiles)))
            log.flush(); print(json.dumps(dict(TRAIN_windows=w['id'], query_profiles=profile_count)), flush=True)
    if profile_count != 432 or len(results) != 27:
        raise ValueError('Incomplete unchanged source readout')
    with (OUT/'references.jsonl').open('x') as log:
        for row in results:
            log.write(json.dumps(row, allow_nan=False)+'\n')
    R.write_summary(OUT, dict(experiment='M28-A16-LEGAL-SUPPORT-DIAGNOSTIC', completed=True,
        protocol=PROTOCOL, source=str(SOURCE.relative_to(R.ROOT)), generator_checkpoint=C.GENERATOR,
        windows=8, reference_window_pairs=27, query_profiles=profile_count, exon_decisions=event_count,
        all_original_geometry_ranks_carry_and_raw_scores_replayed=True, raw_score_max_absolute_difference=maximum_delta,
        reference_balanced_pooled=aggregate(results),
        per_species={s: aggregate([r for r in results if r['species'] == s]) for s in R.ALLOW},
        by_exon_class={k: aggregate([r for r in results if r['exon_class'] == k]) for k in ('one', 'two_to_four', 'five_or_more')},
        by_original_natural_beam_hit={str(hit): aggregate([r for r in results if r['beam_free_hit'] == hit]) for hit in (False, True)},
        reference_conditioned=True, natural_predictions_reused_without_change=True, new_candidate_search=False,
        optimizer_steps_this_run=0, checkpoint_written=False, new_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Repeated TRAIN pressure cohort, forced reference prefixes. Conditional values are diagnostic identities, not a new natural decoder, training result, accuracy, causal proof or novelty.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
