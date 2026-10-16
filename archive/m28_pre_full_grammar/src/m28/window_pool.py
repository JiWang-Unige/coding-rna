"""Fixed reference-free all-window pooling; no change to legacy ownership."""
from .c0_decode import genomic_chain
from .core import select_nonoverlapping


def assemble_window_pool(records, windows, strand):
    """Max existing gain for each complete global chain, then original selector.

    This does not predict new chains or recalibrate scores. The max over windows
    can favor false positives and is not a biological probability.
    """
    windows = [tuple(w) for w in windows]
    if not windows or windows != sorted(set(windows)) or strand not in ('+','-'):
        raise ValueError('Expected nonempty unique sorted windows and a single strand')
    pool = {}; total = 0
    for record in records:
        i = record['window_index']; a,b = windows[i]
        if len(record['chains']) != len(record['score']):
            raise ValueError('One original gain required per candidate')
        for chain,score in zip(record['chains'],record['score']):
            if not chain or not 0 <= chain[0][0] < chain[-1][1] <= b-a:
                raise ValueError('Candidate lies outside its oriented window')
            key = tuple(genomic_chain(chain,a,b,strand)); value = float(score)
            total += 1; old = pool.get(key)
            if old is None or value > old['score'] or (value == old['score'] and i < old['window_index']):
                pool[key] = {'chain':key,'score':value,'window_index':i}
    ordered = [pool[key] for key in sorted(pool)]
    chosen = select_nonoverlapping([r['chain'] for r in ordered],[r['score'] for r in ordered])
    predictions = [ordered[i] for i in chosen]
    return {'predictions':predictions,
            'stages':{'free_unique':sorted(pool),'window_pool_unique':sorted(pool),
                      'selected_unique':[r['chain'] for r in predictions]},
            'counts':{'free_records':total,'free_unique':len(pool),'window_pool_unique':len(pool),
                      'duplicate_free_records':total-len(pool),'selected_unique':len(predictions),
                      'score_field':'score','pool_policy':'max_gain_across_all_generating_windows'}}
