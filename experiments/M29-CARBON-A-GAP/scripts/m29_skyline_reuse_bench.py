"""Timing + identity check: original decode_full_grammar_skyline vs arena-reuse copy on a 1 Mb real strand with Carbon-A offset logits."""
import json, sys, time
import numpy as np, torch
import m29_common as M
from m29_gate_parity import load_tracks, build_logits
from src.m28.full_grammar_skyline import decode_full_grammar_skyline as orig
from m29_skyline_reuse import decode_full_grammar_skyline as reuse
seq = sys.argv[1]; a, n = 6_000_000, int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
rid, sequence, pos, neg = load_tracks(seq)
tt = np.load(M.CACHE / f'tt_{seq}_+.npy')
s = sequence[a:a + n]
lg = torch.from_numpy(build_logits(s, pos[a:a + n], tt[a:a + n], 'A_full'))
out = {}
for _, fn in (('o', orig), ('r', reuse)):
    fn(s[:60000], lg[:60000], 128)      # JIT warm-up (cache=False in the original)
for name, fn in (('original', orig), ('arena_reuse', reuse)):
    t = time.perf_counter(); r = fn(s, lg, 128); dt = time.perf_counter() - t
    out[name] = dict(seconds=dt, s_per_Mb=dt / (n / 1e6), chains=len(r['chains']), gain=r['gain'])
    out[name]['_chains'] = r['chains']
    print(name, {k: v for k, v in out[name].items() if k != '_chains'}, flush=True)
out['identical_chains'] = out['original']['_chains'] == out['arena_reuse']['_chains']
for k in ('original', 'arena_reuse'):
    del out[k]['_chains']
print(out['identical_chains'])
(M.M29 / 'G0b' / f'skyline_arena_bench_{seq}.json').write_text(json.dumps(out, indent=1))
