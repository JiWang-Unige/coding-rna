"""Current-model original-search first-loss collection and qualified gradients.

Only assigned positive queries receive this extra term. Original teacher
supervision remains the trainer's responsibility, including raw-grammar
unsupported references. There is no optimizer, fitting protocol or caller here.
"""
import torch
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.query_prefix_value import START, NEXT_START
from src.m28.query_prefix_value_training import compatible
from trace_frozen_query_beam import trace_query
from analyze_search_cutoff_labels_train import classify, node_at_loss
from query_search_score_replay import replay_node_score, qualified_retention_margin


def raw_support_reason(window, sequence, chain):
    """Report the existing raw/teacher support difference without dropping supervision."""
    sequence = sequence.upper()
    window.validate_chain(chain)
    if len(sequence) != window.length:
        raise ValueError('Original oriented DNA and valid window must agree')
    if len(chain) > 128:
        return 'exon_cap'
    orf, events = ORFIndex(sequence), _end_events(sequence)
    carry, previous_end = '', -1
    for depth, (a,b) in enumerate(chain):
        lower = 0 if depth == 0 else previous_end+Budget().min_intron_bp
        if depth == 0:
            if sequence[a:a+3] != 'ATG':
                return 'start'
        elif a < lower or sequence[a-2:a] != 'AG':
            return 'next_start'
        kind = 'stop' if depth == len(chain)-1 else 'donor'
        choices = _legal_ends(orf,events,a,carry,first_exon=depth == 0)
        match = [c for x,k,c in choices if (x,k) == (b,kind)]
        if len(match) != 1:
            return 'end_'+kind
        carry, previous_end = match[0], b
    return None


def loss_nodes(loss):
    """Represent both same-history local pruning and cross-history pool pruning."""
    if loss['stage'] in ('local_start_topk','local_next_start_topk'):
        target = loss['target']
        chain = tuple(tuple(e) for e in target['chain'])
        event = NEXT_START if chain else START
        return (dict(chain=chain,entry=target['entry'],event=event),
                dict(chain=chain,entry=loss['cutoff_coordinate'],event=event))
    return node_at_loss(loss,True),node_at_loss(loss)


def search_retention_loss(window, sequence, assigned_chains, targets, margin):
    """Observe current natural search, then replay only safe first-loss pairs.

One query per original teacher assignment; no reference path is injected or
kept alive. Loss is averaged over all assigned positive queries, including
returned, protected and unsupported zero terms, not just active violations.
"""
    primary = {tuple(tuple(e) for e in c) for c in targets['positive_chains']}
    zero = window.initial.sum()*0+window.query_logits.sum()*0
    terms, records = [],[]
    for q, reference in sorted(assigned_chains.items()):
        chain = tuple(tuple(e) for e in reference)
        if chain not in primary:
            raise ValueError('Assigned target must remain an original eligible primary')
        reason = raw_support_reason(window,sequence,chain)
        if reason is not None:
            terms.append(zero)
            records.append(dict(query_index=q,chain=chain,raw_grammar_supported=False,
                                unsupported_reason=reason,supervision_disposition='original_teacher_only'))
            continue
        _, observed = trace_query(window,sequence,q,{0:chain})
        trace = observed[0]
        if trace['returned']:
            terms.append(zero)
            records.append(dict(query_index=q,chain=chain,raw_grammar_supported=True,
                                returned=True,first_loss=None,supervision_disposition='already_returned'))
            continue
        target, cutoff = loss_nodes(trace['first_loss'])
        if not compatible(target,chain):
            raise ValueError('Observed target must still be its original positive prefix')
        annotation = classify(cutoff,targets,chain)
        target_score = replay_node_score(window,sequence,q,target)
        cutoff_score = replay_node_score(window,sequence,q,cutoff)
        native_target = trace['first_loss']['target']['path_log_score']
        if abs(float(target_score.detach())-native_target) > 1e-3:
            raise ValueError('Live-graph target replay differs from current original search')
        if trace['first_loss']['stage'] not in ('local_start_topk','local_next_start_topk'):
            key = 'completion_cutoff' if trace['first_loss']['stage'] == 'retained_beam_score_bound' else 'cutoff'
            if abs(float(cutoff_score.detach())-trace['first_loss'][key]['path_log_score']) > 1e-3:
                raise ValueError('Live-graph cutoff replay differs from current original search')
        term = qualified_retention_margin(target_score,cutoff_score,annotation['label'],margin)
        terms.append(term)
        records.append(dict(query_index=q,chain=chain,raw_grammar_supported=True,returned=False,
            first_loss=trace['first_loss'],target_node=target,cutoff_node=cutoff,
            cutoff_annotation=annotation,target_score=float(target_score.detach()),
            cutoff_score=float(cutoff_score.detach()),margin_loss=float(term.detach()),
            supervision_disposition='known_nonmatch_margin' if annotation['label'] == 0 else 'primary_or_unresolved_zero_margin'))
    loss = torch.stack(terms).mean() if terms else zero
    return dict(loss=loss,records=records,positive_queries=len(assigned_chains),
        reference_used_for_observation_only=True,reference_injected_candidates=0,
        beam_width=4,completed_limit_per_query=4,max_exons=128)
