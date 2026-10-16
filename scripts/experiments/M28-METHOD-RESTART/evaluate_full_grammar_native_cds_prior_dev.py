#!/usr/bin/env python3
"""Prepared native-prior DEV evaluation, reusing the unchanged fixed ruler."""
import infer_full_grammar_native_cds_prior_dev as I
import evaluate_full_grammar_dev_baselines as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-COMPARISON'

if __name__ == '__main__':
    I.require_train_closure()
    E.main(I.OUT, I.CHECKPOINT, I.A.IDENTITY, I.A.PROTOCOL, I.METHOD, OUT)
