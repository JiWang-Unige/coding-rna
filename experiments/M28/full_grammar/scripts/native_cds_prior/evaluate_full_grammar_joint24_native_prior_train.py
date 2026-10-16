#!/usr/bin/env python3
"""Unchanged six-TRAIN ruler with the combined-cell paths and identity only."""
import annotate_full_grammar_joint24_native_cds_prior as A
import infer_full_grammar_joint24_native_prior_train as I
import evaluate_full_grammar_original_train as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-DIAGNOSTIC'

if __name__ == '__main__':
    E.main(I.OUT,I.CHECKPOINT,A.IDENTITY,OUT)
