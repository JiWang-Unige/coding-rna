#!/usr/bin/env python3
"""Read-only last-A6-head diagnostic on the original 1491 TRAIN windows."""
import json
import math
import time
from collections import Counter, defaultdict
import numpy as np
import torch
import run_chain_quality as C
from src.m28.labels import build_targets, candidate_targets
from src.m28.training import inject_training

ROOT = C.ROOT
OUT = ROOT/'outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a12_quality_train_diagnostic.md'
COMPONENTS = ('gain', 'additive', 'nonadditive', 'null')


def exon_class(chain):
    return 'one' if len(chain) == 1 else 'two_to_four' if len(chain) <= 4 else 'five_or_more'


def known_mass(labels, weight):
    known = sum(y >= 0 for y in labels)
    return {str(y): weight*labels.count(y)/known if known else 0. for y in (0, 1)}


def overlapping_max(chain, candidates):
    a, b = chain[0][0], chain[-1][1]
    values = [r['gain'] for r in candidates if r['label'] == 0 and
              r['chain'][0][0] < b and r['chain'][-1][1] > a]
    return max(values) if values else None


def profile(rows):
    result = {'candidates': len(rows), 'positive_gain': sum(r['gain'] > 0 for r in rows),
              'nonpositive_gain': sum(r['gain'] <= 0 for r in rows)}
    for key in COMPONENTS:
        values = [r[key] for r in rows]
        result[key] = ({'mean': float(np.mean(values)),
                        'quantiles_0_25_50_75_100': np.quantile(values, [0, .25, .5, .75, 1]).tolist()}
                       if values else None)
    positives = [r for r in rows if r['label'] == 1]
    competing = [r for r in positives if r['max_overlapping_known_nonmatch_gain'] is not None]
    result['positive_candidates_with_known_overlap'] = len(competing)
    result['positive_candidates_not_above_max_known_overlap'] = sum(
        r['gain'] <= r['max_overlapping_known_nonmatch_gain'] for r in competing)
    return result


@torch.inference_mode()
def main():
    fitted = json.loads((C.FIT/'summary.json').read_text())
    C.validate_head(fitted)
    if (C.FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal immutable A6 last head required')
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    C.validate_head(state, state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False)
    head.load_state_dict(state['head'], strict=True)
    del state
    generator = C.load_generator()
    draws, index, windows, catalogs = C.R.Q.load_train_inputs()
    free_rows = C.R.rows(C.FIT/'free_pool.jsonl')
    pools = {r['window']['id']: r for r in free_rows}
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    if (len(free_rows) != 1491 or len(pools) != 1491 or set(pools) != set(windows) or
        sum(len(r['chains']) for r in free_rows) != fitted['pool_counts']['free_candidates']):
        raise ValueError('Complete original free pool must be reused, not regenerated')
    for ident in ids:
        free, w = pools[ident], windows[ident]
        if (free['reference_used_in_pool_generation'] or
            any(w[k] != free['window'][k] for k in C.R.GEOMETRY) or
            w['split'] != 'train' or w['seqid'] not in C.R.ALLOW[w['species']]):
            raise ValueError('Original allowed TRAIN geometry and natural pool required')
    OUT.mkdir()
    began = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
    candidates, window_records = [], {}
    with (OUT/'windows.jsonl').open('x') as log:
        for ordinal, ident in enumerate(ids):
            w, free = windows[ident], pools[ident]
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = C.R.Q.window_inputs(cached, w, 'cuda')
            hidden = generator(features, dna, [w['valid_bases']])['features'][0]
            targets = build_targets(w, catalogs[w['species']], cached['sequence'])
            pool, origins = inject_training(free['chains'], targets['positive_chains'])
            labels = candidate_targets(pool, targets)
            score = head(hidden, pool, w['valid_bases'])
            values = {k: score[k].float().cpu().tolist() for k in ('additive', 'nonadditive', 'logits')}
            null = float(score['null'])
            rows = []
            for chain, origin, label, add, non, gain in zip(
                    pool, origins, labels, values['additive'], values['nonadditive'], values['logits']):
                if (not all(math.isfinite(v) for v in (add, non, gain, null)) or
                    abs(gain-(add+non-null)) > 1e-4 or (origin == 'reference_injected' and label != 1)):
                    raise ValueError('Finite exact component score and unchanged injected positive required')
                rows.append(dict(window_id=ident, species=w['species'], chain=chain, origin=origin,
                                 label=label, exon_class=exon_class(chain), gain=gain,
                                 additive=add, nonadditive=non, null=null))
            for r in rows:
                r['max_overlapping_known_nonmatch_gain'] = overlapping_max(r['chain'], rows) if r['label'] == 1 else None
            known = [r for r in rows if r['label'] >= 0]
            derivative = {str(y): sum(float(torch.sigmoid(torch.tensor(r['gain'])))-y
                                      for r in known if r['label'] == y)/len(known)
                          if known else 0. for y in (0, 1)}
            counts = {'windows': 1, 'free_candidates': len(free['chains']),
                      'reference_injected_candidates': origins.count('reference_injected'),
                      'zero_known_windows': int(not known),
                      **{'label_'+str(y): labels.count(y) for y in (1, 0, -1)}}
            record = dict(window=w, counts=counts, candidates=rows,
                          known_candidate_mean_class_mass=known_mass(labels, 1.),
                          mean_BCE_gain_derivative_by_class=derivative)
            window_records[ident] = record
            candidates.extend(rows)
            log.write(json.dumps(record, allow_nan=False)+'\n')
            if (ordinal+1) % 64 == 0:
                log.flush()
                print(json.dumps({'scored_TRAIN_windows': ordinal+1, 'total': 1491}), flush=True)
            del cached, features, dna, hidden, targets, pool, origins, labels, score
    replay = Counter()
    weighted_mass, weighted_derivative = Counter(), Counter()
    for d in draws:
        r = window_records[d['window_id']]
        replay.update(r['counts'])
        weight = d['importance_weight_uniform_within_stratum']
        weighted_mass.update({y: weight*v for y, v in r['known_candidate_mean_class_mass'].items()})
        weighted_derivative.update({y: weight*v for y, v in r['mean_BCE_gain_derivative_by_class'].items()})
    if any(dict(replay) != e['counts'] for e in fitted['epochs']):
        raise ValueError('Reconstructed draw counts must replay each unchanged A6 epoch')
    groups = defaultdict(list)
    for r in candidates:
        for key in ('all', 'label:'+str(r['label']), 'origin:'+r['origin'],
                    'species:'+r['species']+'|label:'+str(r['label']),
                    'exon:'+r['exon_class']+'|label:'+str(r['label']),
                    'origin:'+r['origin']+'|label:'+str(r['label']),
                    'species:'+r['species']+'|exon:'+r['exon_class']+'|label:'+str(r['label'])):
            groups[key].append(r)
    C.R.write_summary(OUT, dict(experiment='M28-A12-QUALITY-TRAIN-DIAGNOSTIC', completed=True,
        protocol=PROTOCOL, source_generator=C.GENERATOR, source_head=fitted['final_checkpoint'],
        windows=1491, original_draws=1536, unique_candidate_window_pairs=len(candidates),
        replayed_original_draw_counts=dict(replay), profiles={k: profile(v) for k, v in groups.items()},
        original_draw_weighted_class_mass=dict(weighted_mass),
        original_draw_weighted_mean_BCE_gain_derivative_by_class=dict(weighted_derivative),
        reference_used_for_original_free_generation=False, new_candidate_search=False,
        cached_TRAIN_features_used=True, new_backbone_forward=False, DEV_forward_used=False,
        optimizer_steps_this_run=0, head_weights_updated=False, checkpoints_written=False,
        threshold=0., threshold_scan=False, test_setaria_used=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='In-sample last-head TRAIN diagnosis, not accuracy/generalization or causal proof. '
                    'Counts are candidate/window geometry pairs, not independent genes. Label0 is '
                    'current primary-reference-relative nonmatch, not biological truth. No alternative model fitted.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
