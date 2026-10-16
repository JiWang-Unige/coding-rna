#!/usr/bin/env python3
"""Prepared fixed DEV ruler, with strict role-only or actual codon identity."""
import argparse
import torch
import infer_full_grammar_joint_controls_dev as I
import evaluate_full_grammar_dev_baselines as E


def check_candidate(branch, outer, raw, task, species, seqid, length, *, codon_scores=None,
                    checkpoint, identity, protocol):
    train, _ = I.candidate(branch)
    if (not outer['inference_complete']
            or (outer['task'],outer['species'],outer['seqid'],outer['length']) != (task,species,seqid,length)
            or outer['branch'] != branch or outer['checkpoint'] != str(train.CHECKPOINT)
            or checkpoint != train.CHECKPOINT or identity != train.A.IDENTITY or protocol != train.A.PROTOCOL
            or outer['source_fit_job'] != I.SPECS[branch][1] or outer['plan'] != I.PLAN
            or not outer['all_three_TRAIN_flows_closed'] or outer['raw_scope'] != I.I.RAW_SCOPE
            or not outer['DEV_used'] or any(outer[k] for k in
                ('sealed_data_used','reference_used','labels_or_feature_cache_used'))
            or outer['optimizer_steps_this_run'] != 0):
        raise ValueError('Fixed branch, final weight or DEV raw receipt differs')
    if branch == 'CODON_ENERGY':
        if raw['protocol'] != protocol or codon_scores is None:
            raise ValueError('Actual strict final codon parameters and protocol required')
        I.D.validate_raw_report(raw, seqid, length, codon_scores)
        if any(scope['tie_policy'] !=
                'background_continuation_then_state_table_arrival_lower_count_dominates_equal_score'
                for scope in raw['scopes']):
            raise ValueError('Fixed codon tie policy differs')
    else:
        E.check_candidate(outer, raw, task, species, seqid, length,
            checkpoint=checkpoint, identity=identity, protocol=protocol)


def main(branch):
    train, method = I.candidate(branch)
    I.require_train_closure()
    scores = None
    if branch == 'CODON_ENERGY':
        state = torch.load(train.CHECKPOINT, map_location='cpu', weights_only=True)
        train.A.validate_checkpoint(state)
        scores = state['raw_codon_scores']
    def validate(outer, raw, task, species, seqid, length, **expected):
        check_candidate(branch, outer, raw, task, species, seqid, length,
            codon_scores=scores, **expected)
    E.main(I.OUT/branch, train.CHECKPOINT, train.A.IDENTITY, train.A.PROTOCOL,
        method, I.COMPARISON/branch, candidate_validator=validate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch', choices=I.BRANCHES, required=True)
    main(parser.parse_args().branch)
