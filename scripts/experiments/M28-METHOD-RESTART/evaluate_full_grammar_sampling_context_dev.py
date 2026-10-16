#!/usr/bin/env python3
"""Prepared new-candidate comparison, reusing the original frozen DEV ruler."""
import argparse
import infer_full_grammar_sampling_context_dev as I
import evaluate_full_grammar_dev_baselines as E


def main(branch):
    train,method = I.candidate(branch)
    I.require_train_closure()
    E.main(I.OUT/branch,train.CHECKPOINT,train.A.IDENTITY,train.A.PROTOCOL,
        method,I.COMPARISON/branch)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch',choices=I.BRANCHES,required=True)
    main(parser.parse_args().branch)
