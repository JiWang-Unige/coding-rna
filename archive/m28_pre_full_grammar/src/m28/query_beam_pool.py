"""Reference-free completed greedy/beam union for the A6 chain scorer."""


def completed_union(greedy_queries, beam_queries):
    """Keep old completed chains, deduplicate geometry, retain both origins."""
    sources = {}
    for record in greedy_queries:
        if record['terminated']:
            key = tuple(tuple(exon) for exon in record['chain'])
            sources.setdefault(key, set()).add('greedy')
    for query in beam_queries:
        for record in query['completed']:
            key = tuple(tuple(exon) for exon in record['chain'])
            sources.setdefault(key, set()).add('beam4')
    ordered = sorted(sources)
    return [list(chain) for chain in ordered], [sorted(sources[chain]) for chain in ordered]
