"""Internal-validation score (held-out TRAIN chromosomes only, never DEV): pooled exact-chain F1 (CDS>=96) of a run vs the gate on the same chromosomes.
Usage: a2_valscore.py <run_dir> [label]"""
import json, sys
import a2_common as C2
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E
refs = E.Refs()
run = sys.argv[1]
def pooled(r):
    ch = {s: E.load_chains(r, s) for s in C2.HELDOUT}
    res = E.summarize(ch, refs, list(C2.HELDOUT))
    return res
a = pooled(run); g = pooled(K.A1 / 'gate_ownpath')
out = dict(run=run, heldout=list(C2.HELDOUT), run_pooled=a['pooled'], gate_pooled=g['pooled'],
           run_per_contig={s: round(d['F1'], 3) for s, d in a['per_contig'].items()}, gate_per_contig={s: round(d['F1'], 3) for s, d in g['per_contig'].items()},
           delta_F1=a['pooled']['F1'] - g['pooled']['F1'], run_strata={k: round(v['F1'], 2) for k, v in a['by_strata'].items()}, gate_strata={k: round(v['F1'], 2) for k, v in g['by_strata'].items()})
json.dump(out, open(f'{run}/valscore.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('run_pooled', 'gate_pooled', 'delta_F1', 'run_per_contig', 'gate_per_contig')}))
