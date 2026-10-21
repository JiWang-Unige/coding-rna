"""Reproducibility of the re-run: M33 val/dev chains vs the M31 originals (exact-chain F1 and identical chain sets).  Usage: m33_repro_check.py"""
import json
import m33_common as C
import a1_eval as E
refs = E.Refs()
M31 = next(C.K.OUTS.glob('M31*'))
for name, old, seqs in (('val', M31 / 'runB' / 'val_1536', C.HELDOUT), ('dev', M31 / 'devB_final', C.DEV)):
    new = C.OUT / name
    cn = {s: E.load_chains(new, s) for s in seqs}
    co = {s: E.load_chains(old, s) for s in seqs}
    f_new = E.summarize(cn, refs, list(seqs))['pooled']; f_old = E.summarize(co, refs, list(seqs))['pooled']
    same = {s: (len(set(cn[s]) & set(co[s])), len(cn[s]), len(co[s])) for s in seqs}
    print(name, 'F1 new', round(f_new['F1'], 3), 'old', round(f_old['F1'], 3), 'chains common/new/old', same, flush=True)
