"""Seed / pseudo-background inputs of the M36 arms (protocol 22 §2) and the labeler built from them.  No model, no GPU, no reference annotation.
Files:  S2 seeds   DATA/seeds/seeds_S2_<species>.json                (m36_seed_filter.py, CPU preparation)
        S1 seeds   OUT/seeds/seeds_S1_<species>.json                 (m36_s1_seeds.py build; chains with log P >= -1.400 of the zero-shot A run)
        background OUT/seeds/background_<species>.json               ({"contig|strand": [[a, b], ...]}, per-base P(no gene) >= 0.99; used by S1 and S3 only)
Arms: S1 = S1 seeds + background | S2 = S2 seeds, no background | S3 = S2 U S1 (S2 priority) + S1's background | F1 = S1 seeds, every other position background.
"""
import hashlib
import json
from pathlib import Path

import m36_common as W
import m36_seed_labels as SL


def s1_path(species):
    return W.OUT / 'seeds' / f'seeds_S1_{species}.json'


def s2_path(species):
    return W.DATA / 'seeds' / f'seeds_S2_{species}.json'


def background_path(species):
    return W.OUT / 'seeds' / f'background_{species}.json'


def ready_marker(species):
    return W.OUT / 'seeds' / f'READY_{species}'


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def load_background(path):
    raw = json.load(open(path))
    out = {}
    for k, iv in raw.items():
        c, s = k.rsplit('|', 1)
        out[c, s] = [(int(a), int(b)) for a, b in iv]
    return out


def load_arm_inputs(arm, species, s1_file=None, s2_file=None, bg_file=None, contigs=None):
    """-> dict(rows, s1, s2, background, dropped, files {role: (path, md5)}); `rows` are the seeds of the arm (all A contigs; the labeler keeps the training ones)."""
    if arm not in W.ARMS:
        raise ValueError(arm)
    if contigs is None:
        contigs = set(W.TARGETS[species]['A']) if species in W.TARGETS else None
    s1 = s2 = background = None
    files = {}
    if arm in ('S1', 'S3', 'F1'):
        p = Path(s1_file or s1_path(species))
        s1 = SL.load_seed_json(p, contigs)
        files['S1'] = (str(p), md5(p))
    if arm in ('S2', 'S3'):
        p = Path(s2_file or s2_path(species))
        s2 = SL.load_seed_json(p, contigs)
        files['S2'] = (str(p), md5(p))
    if arm in ('S1', 'S3'):                                          # S2 has no pseudo-background (protocol 22 §2), F1 labels every other position background
        p = Path(bg_file or background_path(species))
        background = load_background(p)
        files['background'] = (str(p), md5(p))
    rows, dropped = SL.arm_seeds(arm, s1, s2)
    return dict(arm=arm, species=species, rows=rows, s1=s1, s2=s2, background=background, dropped_s1_overlapping_s2=dropped, files=files)


def make_labeler(inputs, genome_fa=None, lengths=None):
    """SeedLabeler over the species' label-set A TRAINING contigs (validation pseudo-contig and every evaluation-set B contig excluded)."""
    species = inputs['species']
    if lengths is None:
        lens = W.a_lengths(species)
        lengths = {c: lens[c] for c in W.TARGETS[species]['train']}
    return SL.SeedLabeler(species, genome_fa or (W.M1 / species / 'genome.fa'), lengths, inputs['rows'], inputs['background'])


def train_stats(inputs, lab):
    """counts recorded in m36_arm.json"""
    per = {}
    for r in lab.catalog.values():
        per[r['seqid']] = per.get(r['seqid'], 0) + 1
    bg_bp = sum(b - a for iv in (inputs['background'] or {}).values() for a, b in iv)
    return dict(seeds_total=len(inputs['rows']), seeds_train_contigs=len(lab.catalog), seeds_per_contig=per,
                background_intervals=sum(len(v) for v in (inputs['background'] or {}).values()), background_bp_all_A_contigs=bg_bp,
                dropped_s1_overlapping_s2=inputs['dropped_s1_overlapping_s2'])
