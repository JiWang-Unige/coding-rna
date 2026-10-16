"""Finite closed natural competitors for a separate emission-training proposal.

This module has no input loader, search, tensor or training side effects.
The original gold geometries and all original labels remain unchanged.
"""


def augment_training_pool(original_chains, original_labels, natural_chains, natural_labels):
    """Append every new closed natural geometry in its original order.

    Shared geometries require identical frozen labels. New positives would
    change the original supervision domain and are rejected, not injected.
    Unknown additions remain -1; the existing structured loss excludes them,
    whereas subsequent reference-free prediction must not filter them out.
    """
    if (len(original_chains) != len(original_labels) or
            len(natural_chains) != len(natural_labels)):
        raise ValueError('One saved label per geometry required')
    if any(y not in (-1, 0, 1) for y in (*original_labels, *natural_labels)):
        raise ValueError('Original 1/0/-1 labels required')
    original_keys = [tuple(map(tuple, chain)) for chain in original_chains]
    natural_keys = [tuple(map(tuple, chain)) for chain in natural_chains]
    if (len(set(original_keys)) != len(original_keys) or
            len(set(natural_keys)) != len(natural_keys)):
        raise ValueError('Closed original and natural geometries must be unique')
    by_geometry = dict(zip(original_keys, original_labels))
    chains = [[list(exon) for exon in chain] for chain in original_chains]
    labels = list(original_labels)
    added = []
    for i, (geometry, label) in enumerate(zip(natural_keys, natural_labels)):
        if geometry in by_geometry:
            if by_geometry[geometry] != label:
                raise ValueError('Common geometry frozen labels differ')
            continue
        if label == 1:
            raise ValueError('Natural positive absent from the original gold domain')
        chains.append([list(exon) for exon in geometry])
        labels.append(label)
        added.append(i)
    return dict(chains=chains, labels=labels, added_natural_indices=added)
