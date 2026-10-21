import json, sys, numpy as np, torch
import a1_common as K, a1_model as MD, a1_offset as O, a1_dp as DP, a1_train as T
import prepare_full_grammar_fresh_draws as PP, fit_full_grammar_objective_contrast as C
plan = [json.loads(x) for x in (K.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
windows, strata, catalogs, original = PP.load_population()
sequences, _ = C.P.read_train_DNA()
tracks = {s: O.load_tracks(s) for s in {w['seqid'] for w in windows.values()}}
DP._install()
for i in range(int(sys.argv[1]), int(sys.argv[2])):
    r = plan[i]; w = windows[r['window_id']]
    chromosome, geometry, context, offset, targets, counts = T.example(w, sequences, catalogs)
    a, b = geometry['genomic_context']
    off = O.context_offset(tracks[w['seqid']], w['strand'], a, b, len(chromosome))
    seg = np.asarray(targets['seg'])[:w['valid_bases']]
    res = DP.window_loss_grad(context, off, offset, offset + w['valid_bases'], {k: np.asarray(targets[k]) for k in ('seg','region')})
    print(i, w['species'][:3], w['strand'], 'fine', counts[0], 'coarse', counts[1], 'gene_roles', int((seg >= 1).sum()), 'bg', int((seg == 0).sum()),
          'loss', res['loss'], 'logZall', res['log_Z_all'], 'logZobs', res['log_Z_obs'], 'gmax', float(np.abs(res['grad']).max()), 'dp_s', round(res['seconds'], 1), flush=True)
