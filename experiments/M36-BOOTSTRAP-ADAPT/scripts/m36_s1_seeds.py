"""S1 seeds + pseudo-background from the zero-shot run on the label-set A contigs (protocol 22 §2, rows S1 / S3 / F1).  CPU only; no reference annotation is read.
Inputs (per A contig and strand, written by m36_cpu_queue.py from the zero-shot run):
  post_A/tau1/<contig>_<strand>.json        m33_posterior.py: log P(g in gene set) of every decoded chain at tau = 1 (M33 frozen selection: tau* = 1, thr* = -1.400)
  nogene_A/tau1/nogene_<contig>_<strand>.npy (+ .json)   m36_nogene_marginal.py: per-base P(no gene covers the base) (oriented float32)
Outputs (OUT/seeds/):
  seeds_S1_<species>.json        chains with log P >= -1.400270636 (P >= 0.2466), CDS >= 96 bp (the project ruler's chains), passing the frozen chain_checks; seed-row format
  seeds_S1_P90_<species>.json    sensitivity set P >= 0.9 (not used by any arm)
  background_<species>.json      {"contig|strand": [[a, b], ...]} genomic 0-based half-open runs with P(no gene) >= 0.99 (`pseudo_background(..., 0.99)`)
  seed_stats_S1_<species>.json   funnel / counts per contig / exon-count strata / background coverage;   READY_<species> marker (written last)
usage: m36_s1_seeds.py build <species> [--post DIR] [--nogene DIR] [--out DIR]
Not a model-dependent choice: the threshold, tau and the 96-bp rule are the M33 / project constants in m36_common.
"""
import argparse
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

import m36_common as W
import m36_seed_labels as SL
import m36_train_inputs as TI
import manifest as MF


def chains_from_posterior(post, thr=W.S1_LOGP_THR, min_cds_bp=W.S1_MIN_CDS_BP):
    """one posterior JSON ({seq, strand, genes: [{exons, gain, log_p}]}) -> (kept [(exons, log_p)], funnel counts).  exons are genomic, ascending."""
    f = Counter(decoded=len(post['genes']))
    keep = []
    for g in post['genes']:
        ex = sorted((int(a), int(b)) for a, b in g['exons'])
        lp = float(g['log_p'])
        if lp < thr:
            continue
        f['logp_ge_thr'] += 1
        if sum(b - a for a, b in ex) < min_cds_bp:
            f['dropped_cds_lt_min'] += 1
            continue
        keep.append((ex, lp))
    f['kept_before_checks'] = len(keep)
    return keep, f


def make_rows(species, contig, strand, kept, dna, idx0, label='S1'):
    """kept: [(exons, log_p)] -> seed rows passing the frozen chain_checks, ids {prefix}{label}_{n:06d}; -> (rows, fail Counter)"""
    rows, fails = [], Counter()
    for k, (ex, lp) in enumerate(kept):
        r = SL.seed_row(contig, strand, SL.phases(ex, strand), f'{W.PREFIX.get(species, species[:3])}{label}_{idx0 + len(rows):06d}',
                        log_p=lp, P=math.exp(lp), source='zero-shot decoded chain')
        failures, motifs = MF.chain_checks(r, dna)
        if failures:
            fails.update(failures)
            fails['dropped'] += 1
            continue
        rows.append(r)
    return rows, fails


def build(species, post_dir, nogene_dir, out_dir, thr=W.S1_LOGP_THR, genome=None, lens=None, train=None):
    """lens / train: override the label-set A contig lengths / training contigs (tests with synthetic data only)"""
    t0 = time.time()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lens = lens or W.a_lengths(species)
    contigs = list(lens)
    if species in W.TARGETS:
        W.assert_not_B(species, contigs)
    seqs = MF.read_allowed(Path(genome) if genome else W.M1 / species / 'genome.fa', lens)
    s1, p90, bg_all = [], [], {}
    funnel = Counter()
    per_contig = {}
    bg_stats = {}
    for c in contigs:
        for st in '+-':
            pj = Path(post_dir) / f'tau{W.S1_TAU:g}' / f'{c}_{st}.json'
            post = json.load(open(pj))
            assert post['seq'] == c and post['strand'] == st and abs(post['tau'] - W.S1_TAU) < 1e-12, (pj, post['seq'], post['strand'], post['tau'])
            kept, f = chains_from_posterior(post, thr)
            kept90, _ = chains_from_posterior(post, W.S1_LOGP_P90)
            rows, fails = make_rows(species, c, st, kept, seqs[c], len(s1))
            rows90, fails90 = make_rows(species, c, st, kept90, seqs[c], len(p90), label='S1P90')
            s1 += rows
            p90 += rows90
            f.update(chain_check_failures=sum(fails.values()) - fails['dropped'], chain_check_dropped=fails['dropped'], kept=len(rows))
            funnel.update(f)
            per_contig[f'{c}|{st}'] = dict(decoded=f['decoded'], logp_ge_thr=f['logp_ge_thr'], kept=len(rows), P_ge_0_9=len(rows90), dropped_checks=dict(fails),
                                           log_partition=post['log_partition'])
            # pseudo-background
            npy = Path(nogene_dir) / f'tau{W.S1_TAU:g}' / f'nogene_{c}_{st}.npy'
            p = np.load(npy)
            iv = SL.pseudo_background(p, st, lens[c], W.PSEUDO_BACKGROUND_MIN_P)
            js = json.load(open(npy.with_suffix('.json')))
            assert [list(x) for x in iv] == [list(x) for x in js['pseudo_background_genomic']], f'{npy}: json segments differ from pseudo_background()'
            bg_all[f'{c}|{st}'] = [[int(a), int(b)] for a, b in iv]
            bg_stats[f'{c}|{st}'] = dict(intervals=len(iv), bp=int(sum(b - a for a, b in iv)), fraction=float(sum(b - a for a, b in iv) / lens[c]), mean_P_nogene=float(p.mean()))
    train = set(train if train is not None else W.TARGETS[species]['train'])
    train_bp = sum(lens[c] for c in train)
    n_train = sum(1 for r in s1 if r['seqid'] in train)
    strata = lambda rows: dict(Counter('1' if r['n_exons'] == 1 else '2-5' if r['n_exons'] <= 5 else '>=6' for r in rows))
    # a seed must not sit inside a pseudo-background run of its own strand (P(gene) >= 0.247 vs P(no gene) >= 0.99 are incompatible); count violations (expected 0)
    overlap = 0
    for r in s1:
        for a, b in bg_all.get(f"{r['seqid']}|{r['strand']}", []):
            if a < r['span'][1] and b > r['span'][0]:
                overlap += 1
                break
    stats = dict(species=species, tau=W.S1_TAU, threshold_log_p=thr, threshold_P=math.exp(thr), min_cds_bp=W.S1_MIN_CDS_BP, funnel=dict(funnel),
                 seeds_S1=len(s1), seeds_S1_train_contigs=n_train, seeds_per_Mb_train=1e6 * n_train / train_bp, seeds_S1_P_ge_0_9=len(p90),
                 seeds_S1_P_ge_0_9_train_contigs=sum(1 for r in p90 if r['seqid'] in train), exon_strata_all_A=strata(s1), exon_strata_P90=strata(p90),
                 seeds_overlapping_own_strand_background=overlap, background_total_bp=sum(v['bp'] for v in bg_stats.values()),
                 background_fraction_of_A_strands=sum(v['bp'] for v in bg_stats.values()) / (2 * sum(lens.values())), per_contig_strand=per_contig,
                 background_per_contig_strand=bg_stats, contigs=lens, train_contigs=sorted(train), post_dir=str(post_dir), nogene_dir=str(nogene_dir),
                 reference_used=False, seconds=round(time.time() - t0, 1))
    for name, obj in ((f'seeds_S1_{species}.json', s1), (f'seeds_S1_P90_{species}.json', p90), (f'background_{species}.json', bg_all),
                      (f'seed_stats_S1_{species}.json', stats)):
        tmp = out_dir / (name + '.tmp')
        json.dump(obj, open(tmp, 'w'))
        tmp.replace(out_dir / name)
    (out_dir / f'READY_{species}').write_text(json.dumps(dict(time=time.strftime('%Y-%m-%dT%H:%M:%S'), seeds_S1=len(s1), seeds_S1_train=n_train)))
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('build',))
    ap.add_argument('species', choices=W.SPECIES)
    ap.add_argument('--post', default=str(W.OUT / 'post_A'))
    ap.add_argument('--nogene', default=str(W.OUT / 'nogene_A'))
    ap.add_argument('--out', default=str(W.OUT / 'seeds'))
    a = ap.parse_args()
    W.K.repoint()
    st = build(a.species, a.post, a.nogene, a.out)
    print(json.dumps({k: st[k] for k in ('species', 'seeds_S1', 'seeds_S1_train_contigs', 'seeds_per_Mb_train', 'seeds_S1_P_ge_0_9', 'background_fraction_of_A_strands',
                                         'seeds_overlapping_own_strand_background', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
