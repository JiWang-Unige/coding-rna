"""M31-C1-NTV2 shared data/sampler helpers (frozen pipeline code is imported, never edited).

TRAIN population = the original-six TRAIN chromosomes of M28 (Arabidopsis NC_003070.9 / NC_003071.7 / NC_003076.8, rice
NC_089043.1 / NC_089044.1 / NC_089046.1).  Following M29-A2 the two chromosomes NC_003071.7 and NC_089044.1 are INTERNAL
VALIDATION chromosomes: never sampled for training, never DEV.  DEV (NC_003074.8, NC_089041.1) is used exactly once, at the end.
"""
import math
import os
import re
import sys
from pathlib import Path

import ntv2_adapter as A                                   # sets sys.path (pipeline/m28, M29-A1, M29-CARBON-A-GAP, ...)
import a1_common as K                                      # noqa: F401  (M29-A1 helper: repoint of Baobab-absolute Paths)
import numpy as np
import torch

import fit_full_grammar_objective_contrast as C             # frozen training library
import prepare_full_grammar_fresh_draws as PP               # frozen sampler / population
K.repoint()
G, Q = C.G, C.Q

HELDOUT = ('NC_003071.7', 'NC_089044.1')
SUBSTITUTE = {'NC_003071.7': 'NC_003070.9', 'NC_089044.1': 'NC_089043.1'}   # stratum remap so the frozen schedule keeps its structure
TRAIN_SEQS = ('NC_003070.9', 'NC_003071.7', 'NC_003076.8', 'NC_089043.1', 'NC_089044.1', 'NC_089046.1')
DEV_SEQS = ('NC_003074.8', 'NC_089041.1')
OUT = A.M31


def lr_at(step, total, lr_peak, lr_final, warmup_frac=.05):
    w = max(1, int(math.ceil(warmup_frac * total)))
    if step < w:
        return lr_peak * (step + 1) / w
    t = (step - w) / max(1, total - w)
    return lr_final + .5 * (lr_peak - lr_final) * (1 + math.cos(math.pi * t))


def example(w, sequences, catalogs):
    """Copy of the frozen `C.example` minus the comparison with a stored plan row.  Returns the genomic context geometry
    of the frozen 24,576-bp core +- one 12,288 neighbour, the oriented context string, the core offset and the core targets."""
    chromosome = sequences[w['species']][w['seqid']]
    geometry = G.context_geometry(w['start'], w['end'], len(chromosome))
    a, b = geometry['genomic_context']
    raw = chromosome[w['start']:w['end']]
    core = raw if w['strand'] == '+' else Q.M.rc(raw)
    targets = C.P.P.build_targets(w, catalogs[w['species']], core)
    n = w['valid_bases']
    fine, region = (torch.as_tensor(targets[k])[:n] for k in ('seg', 'region'))
    counts = (int((fine >= 0).sum()), int(((fine < 0) & (region >= 0)).sum()), int(((fine < 0) & (region < 0)).sum()))
    offset = w['start'] - a if w['strand'] == '+' else b - w['end']
    context = re.sub('[^ACGT]', 'N', chromosome[a:b])
    if w['strand'] == '-':
        context = Q.M.rc(context)
    return chromosome, geometry, context, offset, targets, counts


def schedule(total_draws, strata, original, sampling_seed=1):
    """Frozen fresh-draw schedule (seed 1, within-stratum uniform with replacement) with the internal-validation chromosomes
    remapped to a same-species/strand training chromosome stratum -> HELDOUT chromosomes are never sampled."""
    passes = int(math.ceil(total_draws / len(original)))
    strata2 = dict(strata)
    for (species, seqid, strand) in list(strata):
        if seqid in HELDOUT:
            strata2[(species, seqid, strand)] = strata[(species, SUBSTITUTE[seqid], strand)]
    return PP.fresh_schedule(original, strata2, passes=passes, seed=sampling_seed)[:total_draws]


def atomic_save(obj, path):
    tmp = str(path) + '.tmp'
    torch.save(obj, tmp)
    os.replace(tmp, path)


def parse_cpus(spec):
    if not spec:
        return None
    cpus = []
    for part in spec.split(','):
        lo, _, hi = part.partition('-')
        cpus += list(range(int(lo), int(hi or lo) + 1))
    return cpus
