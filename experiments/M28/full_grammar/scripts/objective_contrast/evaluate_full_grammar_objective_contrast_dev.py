#!/usr/bin/env python3
"""Prepared DEV comparison; unchanged native assets, ruler and denominators."""
import argparse
import infer_full_grammar_objective_contrast_dev as I
import evaluate_full_grammar_dev_baselines as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-COMPARISON'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=I.A.ARMS, required=True)
    args = parser.parse_args()
    arm = args.arm
    I.require_train_closure()
    E.main(I.OUT/arm, I.final_checkpoint(arm), dict(I.A.IDENTITY, arm=arm),
           I.A.PROTOCOL, 'M28_'+arm+'_final9216', OUT/arm)
