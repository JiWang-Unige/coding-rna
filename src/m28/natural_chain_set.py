"""A15 natural-domain projection of the closed A12 TRAIN records."""


def natural_record(record):
    source = record['candidates']
    if any(r['origin'] not in ('free', 'reference_injected') for r in source):
        raise ValueError('Expected the closed free/injected candidate domain')
    local = [dict(r) for r in source if r['origin'] == 'free']
    labels = [r['label'] for r in local]
    return {**record, 'candidates': local, 'counts': dict(
        windows=1, free_candidates=len(local), reference_injected_candidates=0,
        zero_known_windows=int(not any(y >= 0 for y in labels)),
        **{'label_'+str(y): labels.count(y) for y in (1, 0, -1)})}
