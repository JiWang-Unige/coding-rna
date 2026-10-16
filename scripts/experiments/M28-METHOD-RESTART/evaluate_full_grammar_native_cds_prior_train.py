#!/usr/bin/env python3
"""Reuse the unchanged TRAIN ruler, with explicit native-prior prediction paths."""
import annotate_full_grammar_native_cds_prior as A
import infer_full_grammar_native_cds_prior_train as I
import evaluate_full_grammar_original_train as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-DIAGNOSTIC'

if __name__ == '__main__':
    E.main(I.OUT, I.CHECKPOINT, A.IDENTITY, OUT)
