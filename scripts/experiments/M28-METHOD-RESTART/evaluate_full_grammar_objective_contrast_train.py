#!/usr/bin/env python3
"""Unchanged original TRAIN ruler, explicit arm output and final checkpoint."""
import argparse
import annotate_full_grammar_objective_contrast as A
import infer_full_grammar_objective_contrast_train as I
import evaluate_full_grammar_original_train as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=A.ARMS, required=True)
    args = parser.parse_args()
    arm = args.arm
    checkpoint = I.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST'/arm/'step_009216.pt'
    E.main(I.OUT/arm, checkpoint, dict(A.IDENTITY, arm=arm), OUT/arm)
