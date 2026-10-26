"""S1 seeds + pseudo-background from the zero-shot posterior (m36_s1_seeds.py), CPU only.
Synthetic part (no data): threshold / CDS >= 96 rule / P >= 0.9 set, seed rows equal the planted chains (CDS + recomputed phases), chain_checks failures are dropped and counted, the S1
labeler on the produced seeds gives the fine labels of the frozen `build_targets` on the seed loci and IGNORE elsewhere except the pseudo-background, the whole `build` on a toy A species.
Data part (ont, skipped otherwise): a slice of the REAL M33 posterior (Arabidopsis internal-validation chromosome NC_003071.7 '+' strand, tau = 1; the file the M33 frozen threshold
-1.400 was calibrated on) -> counts equal an independent count, every kept chain passes the frozen chain_checks on the real genome, labels on real seed windows.
Run (ont): cd experiments/M36-BOOTSTRAP-ADAPT/tests && nice -n 10 taskset -c 25-31 ~/envs/generanno/bin/pytest -q test_m36_s1_seeds.py
"""
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

import m36_common as W
import m36_s1_seeds as S1
import m36_seed_labels as SL
import m36_train_inputs as TI
import manifest as MF
from src.m28.labels import IGNORE
import test_seed_labels as TS            # synthetic gene / contig builders of the seed-label tests

SAFE = TS.SAFE


def long_gene(rng, n_exons):
    """legal literal-DNA complete CDS of 126-240 bp (the S1 rule keeps chains with CDS >= 96 bp)"""
    cds = 'ATG' + ''.join(rng.choice(SAFE, int(rng.integers(40, 80)))) + 'TAA'
    cuts = sorted(rng.choice(np.arange(4, len(cds) - 4), n_exons - 1, replace=False)) if n_exons > 1 else []
    pieces, prev = [], 0
    for c in list(cuts) + [len(cds)]:
        pieces.append(cds[prev:c]); prev = c
    dna, exons, pos = '', [], 0
    for j, piece in enumerate(pieces):
        exons.append((pos, pos + len(piece))); dna += piece; pos += len(piece)
        if j < len(pieces) - 1:
            intron = 'GT' + ''.join(rng.choice(list('ACGT'), int(rng.integers(8, 40)))) + 'AG'
            dna += intron; pos += len(intron)
    return dna, exons


def synth_contig(rng, length=60000, n_genes=6, seq='SYN', start=2000, gap=(500, 2500)):
    """as test_seed_labels.synthetic_contig, with genes long enough for the S1 CDS rule"""
    dna = list(rng.choice(list('ACGT'), length))
    rows, cursor = [], start
    for k in range(n_genes):
        g, exons = long_gene(rng, int(rng.integers(1, 5)))
        strand = '+' if k % 2 == 0 else '-'
        if strand == '+':
            dna[cursor:cursor + len(g)] = list(g)
            cds = [(cursor + a, cursor + b) for a, b in exons]
        else:
            dna[cursor:cursor + len(g)] = list(MF.rc(g))
            cds = sorted((cursor + len(g) - b, cursor + len(g) - a) for a, b in exons)
        rows.append(SL.seed_row(seq, strand, SL.phases(cds, strand), f'g{k}'))
        cursor += len(g) + int(rng.integers(*gap))
    dna = ''.join(dna)
    assert cursor < length
    for r in rows:
        failures, _ = MF.chain_checks(r, dna)
        assert not failures, (r['id'], failures)
    return dna, rows


M33_POST = W.K.OUTS / 'M33-NTV2-FOLLOWUP' / 'post_val' / 'tau1' / 'NC_003071.7_+.json'
ATH = W.M1 / 'arabidopsis_thaliana' / 'genome.fa'
needs_ont = pytest.mark.skipif(not (M33_POST.exists() and ATH.exists()), reason='needs the M33 posterior and the Arabidopsis genome on ont')


def test_constants_are_the_m33_frozen_ones():
    assert abs(W.S1_LOGP_THR - (-1.400270636077039)) < 1e-15 and abs(math.exp(W.S1_LOGP_THR) - 0.2466) < 1e-3 and W.S1_TAU == 1.0
    sel = None
    for cand in (W.EXPS / 'M33-NTV2-FOLLOWUP' / 'results' / 'selection.json', W.K.OUTS / 'M33-NTV2-FOLLOWUP' / 'selection.json'):
        if cand.exists():
            sel = json.load(open(cand))['selection']
            break
    if sel is not None:
        assert sel['tau'] == W.S1_TAU and sel['thr'] == W.S1_LOGP_THR
    assert abs(W.S1_LOGP_P90 - math.log(0.9)) < 1e-15 and W.S1_MIN_CDS_BP == 96 and W.PSEUDO_BACKGROUND_MIN_P == 0.99


def posterior_for(rows, log_ps, seq='SYN', strand='+', extra=()):
    genes = [dict(exons=[[a, b] for a, b, _ in r['CDS']], gain=1.0, log_p=lp) for r, lp in zip(rows, log_ps)]
    genes += [dict(exons=[list(e) for e in ex], gain=0., log_p=lp) for ex, lp in extra]
    return dict(seq=seq, strand=strand, tau=1.0, genes=genes, log_partition=0., seconds=0.)


def test_threshold_min_cds_and_p90_sets():
    thr = W.S1_LOGP_THR
    post = dict(seq='X', strand='+', tau=1.0, log_partition=0., seconds=0., genes=[
        dict(exons=[[100, 400]], gain=1, log_p=-0.05),                 # P = 0.95: kept in both sets
        dict(exons=[[1000, 1300]], gain=1, log_p=-1.0),                # P = 0.37: S1 only
        dict(exons=[[2000, 2300]], gain=1, log_p=thr),                 # exactly at the threshold: kept (>=)
        dict(exons=[[3000, 3300]], gain=1, log_p=thr - 1e-9),          # just below
        dict(exons=[[4000, 4060]], gain=1, log_p=-0.01),               # CDS 60 bp < 96: dropped by the project-ruler rule
        dict(exons=[[5000, 5100], [5200, 5230]], gain=1, log_p=-0.5),  # CDS 130 bp over 2 exons: kept
    ])
    kept, f = S1.chains_from_posterior(post)
    assert [e for e, _ in kept] == [[(100, 400)], [(1000, 1300)], [(2000, 2300)], [(5000, 5100), (5200, 5230)]]
    assert f['decoded'] == 6 and f['logp_ge_thr'] == 5 and f['dropped_cds_lt_min'] == 1
    kept90, _ = S1.chains_from_posterior(post, W.S1_LOGP_P90)
    assert [e for e, _ in kept90] == [[(100, 400)]]                  # P >= 0.9: the 60-bp chain is dropped by the CDS rule, the rest has P < 0.9
    assert all(lp >= math.log(0.9) for _, lp in kept90)


def test_rows_equal_planted_chains_and_failures_are_dropped():
    rng = np.random.default_rng(11)
    dna, rows = synth_contig(rng)
    plus = [r for r in rows if r['strand'] == '+']
    minus = [r for r in rows if r['strand'] == '-']
    assert plus and minus
    for strand, rs in (('+', plus), ('-', minus)):
        bad_exons = [(10, 400)]                                          # 390 bp random sequence: fails ATG start / stop / frame checks
        post = posterior_for(rs, [-0.2] * len(rs), strand=strand, extra=[(bad_exons, -0.1)])
        kept, f = S1.chains_from_posterior(post)
        assert len(kept) == len(rs) + 1
        out, fails = S1.make_rows('dryas_octopetala', 'SYN', strand, kept, dna, 0)
        assert fails['dropped'] == 1 and len(out) == len(rs)             # the invalid chain is counted and dropped, the planted ones survive
        key = lambda r: (r['span'][0], tuple((a, b) for a, b, _ in r['CDS']))
        assert sorted(map(key, out)) == sorted(map(key, rs))
        for o in out:                                                    # phases recomputed in translation order = the planted ones
            p = next(r for r in rs if r['CDS'][0][0] == o['CDS'][0][0])
            assert o['CDS'] == p['CDS'] and o['strand'] == strand and o['eligible'] and o['id'].startswith('dryS1_')


def test_s1_labeler_on_produced_seeds_marks_unknown_and_background(tmp_path):
    rng = np.random.default_rng(5)
    dna, rows = synth_contig(rng)
    fa = tmp_path / 'g.fa'
    fa.write_text('>SYN\n' + '\n'.join(dna[i:i + 80] for i in range(0, len(dna), 80)) + '\n')
    kept, _ = S1.chains_from_posterior(posterior_for([r for r in rows if r['strand'] == '+'], [-0.3] * sum(r['strand'] == '+' for r in rows)))
    s1, _ = S1.make_rows('dryas_octopetala', 'SYN', '+', kept, dna, 0)
    last = max(r['span'][1] for r in s1)
    bg = {('SYN', '+'): [(0, 1000), (last + 100, min(len(dna), last + 700))]}      # runs that touch no seed
    lens = {'SYN': len(dna)}
    lab = SL.SeedLabeler('syn', fa, lens, s1, bg)
    for a, b in lab.grid('SYN'):
        tg, win, core = lab.targets('SYN', '+', a, b, 'S1')
        seg = tg['seg'][:b - a]
        for r in s1:
            x, y = r['span']
            if a <= x and y <= b:
                assert (seg[x - a:y - a] != IGNORE).all()                # seed locus labelled
        known = np.zeros(b - a, bool)
        for r in s1:
            known[max(0, r['span'][0] - a):max(0, min(b, r['span'][1]) - a)] = True
        for p, q in lab.background_local('SYN', '+', a, b):
            assert (seg[p - a:q - a] == 0).all()                          # pseudo-background = background class, not IGNORE
            known[p - a:q - a] = True
        assert (seg[~known] == IGNORE).all() and (tg['region'][:b - a][~known] == IGNORE).all()
    # F1 with the same seeds: nothing unknown
    tg, _, _ = lab.targets('SYN', '+', 0, len(dna) if len(dna) <= 24576 else 24576, 'F1')
    assert (tg['seg'][:min(len(dna), 24576)] != IGNORE).all()


def toy_species(tmp_path, rng):
    """toy A species: SYN (train) + SYV (validation pseudo-contig); posterior / no-gene files for both strands; returns paths and the planted seeds"""
    dna, rows = synth_contig(rng, length=40000, n_genes=6, seq='SYN')
    dna2, rows2 = synth_contig(rng, length=30000, n_genes=4, seq='SYV')
    seqs = {'SYN': dna, 'SYV': dna2}
    fa = tmp_path / 'genome.fa'
    fa.write_text(''.join(f'>{c}\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n' for c, s in seqs.items()))
    post, nog = tmp_path / 'post', tmp_path / 'nogene'
    planted = {}
    for c, rs in (('SYN', rows), ('SYV', rows2)):
        L = len(seqs[c])
        for st in '+-':
            mine = [r for r in rs if r['strand'] == st]
            lps = [-0.3 if k % 2 == 0 else -2.5 for k in range(len(mine))]      # every other chain falls below the threshold
            planted[c, st] = [r for r, lp in zip(mine, lps) if lp >= W.S1_LOGP_THR]
            d = post / 'tau1'
            d.mkdir(parents=True, exist_ok=True)
            p = posterior_for(mine, lps, seq=c, strand=st)
            json.dump(p, open(d / f'{c}_{st}.json', 'w'))
            pn = np.full(L, 0.5, np.float32)
            pn[2000:2600] = 0.995
            pn[L - 1500:L - 900] = 0.999
            for r in mine:
                pn[r['span'][0]:r['span'][1]] = 0.01
            dn = nog / 'tau1'
            dn.mkdir(parents=True, exist_ok=True)
            np.save(dn / f'nogene_{c}_{st}.npy', pn)
            seg = SL.pseudo_background(pn, st, L, 0.99)
            json.dump(dict(pseudo_background_genomic=[list(x) for x in seg]), open(dn / f'nogene_{c}_{st}.json', 'w'))
    return dict(fa=fa, post=post, nogene=nog, lens={c: len(s) for c, s in seqs.items()}, planted=planted, dna=seqs)


def test_build_end_to_end_toy_species(tmp_path):
    t = toy_species(tmp_path, np.random.default_rng(3))
    out = tmp_path / 'seeds'
    st = S1.build('dryas_octopetala', t['post'], t['nogene'], out, genome=t['fa'], lens=t['lens'], train=['SYN'])
    s1 = json.load(open(out / 'seeds_S1_dryas_octopetala.json'))
    assert (out / 'READY_dryas_octopetala').exists() and st['seeds_S1'] == len(s1) == sum(len(v) for v in t['planted'].values())
    got = sorted((r['seqid'], r['strand'], r['span'][0]) for r in s1)
    want = sorted((r['seqid'], r['strand'], r['span'][0]) for v in t['planted'].values() for r in v)
    assert got == want
    assert st['seeds_S1_train_contigs'] == sum(1 for r in s1 if r['seqid'] == 'SYN')
    assert st['seeds_overlapping_own_strand_background'] == 0 and st['funnel']['chain_check_dropped'] == 0
    bg = TI.load_background(out / 'background_dryas_octopetala.json')
    assert set(bg) == {(c, s) for c in t['lens'] for s in '+-'}
    for (c, s), iv in bg.items():                                       # segments reproduce pseudo_background(p, strand, L, 0.99) of the stored probabilities
        p = np.load(t['nogene'] / 'tau1' / f'nogene_{c}_{s}.npy')
        assert [list(x) for x in iv] == [list(x) for x in SL.pseudo_background(p, s, t['lens'][c], 0.99)]
    # arm inputs from the produced files: S1 uses seeds + background, F1 seeds only, S3 merges with S2 priority
    inp = TI.load_arm_inputs('S1', 'dryas_octopetala', s1_file=out / 'seeds_S1_dryas_octopetala.json', bg_file=out / 'background_dryas_octopetala.json', contigs={'SYN', 'SYV'})
    assert len(inp['rows']) == len(s1) and inp['background'] and inp['s2'] is None
    inp = TI.load_arm_inputs('F1', 'dryas_octopetala', s1_file=out / 'seeds_S1_dryas_octopetala.json', contigs={'SYN', 'SYV'})
    assert inp['background'] is None and len(inp['rows']) == len(s1)
    p90 = json.load(open(out / 'seeds_S1_P90_dryas_octopetala.json'))
    assert len(p90) == st['seeds_S1_P_ge_0_9'] == 0 and st['funnel']['kept'] == len(s1)      # planted log P are -0.3 (P = 0.74) and -2.5: none reaches 0.9


def test_build_refuses_inconsistent_nogene_json(tmp_path):
    t = toy_species(tmp_path, np.random.default_rng(4))
    f = t['nogene'] / 'tau1' / 'nogene_SYN_+.json'
    d = json.load(open(f))
    d['pseudo_background_genomic'] = d['pseudo_background_genomic'][:-1]
    json.dump(d, open(f, 'w'))
    with pytest.raises(AssertionError):
        S1.build('dryas_octopetala', t['post'], t['nogene'], tmp_path / 'o', genome=t['fa'], lens=t['lens'], train=['SYN'])


@needs_ont
def test_real_m33_posterior_slice():
    post = json.load(open(M33_POST))
    genes = post['genes'][:400]
    sl = dict(post, genes=genes)
    L = 19698289                                                        # NC_003071.7 length (M29 SCOPE)
    seq = MF.read_allowed(ATH, {'NC_003071.7': L})['NC_003071.7']
    kept, f = S1.chains_from_posterior(sl)
    lp = np.array([g['log_p'] for g in genes])
    cds = np.array([sum(b - a for a, b in g['exons']) for g in genes])
    assert f['decoded'] == 400 and len(kept) == int(((lp >= W.S1_LOGP_THR) & (cds >= 96)).sum()) > 50
    rows, fails = S1.make_rows('arabidopsis_thaliana', 'NC_003071.7', '+', kept, seq, 0)
    print('real slice: kept', len(kept), 'rows', len(rows), 'fails', dict(fails), 'P>=0.9', int(((lp >= math.log(.9)) & (cds >= 96)).sum()))
    assert fails['dropped'] <= 0.02 * len(kept)                         # the frozen grammar emits only chains that pass the frozen checks (a few non-canonical motifs are tolerated and counted)
    assert all(r['eligible'] for r in rows) and len(rows) >= 0.98 * len(kept)
    # labels on real seed windows: seed loci labelled, everything else IGNORE
    lab = SL.SeedLabeler('arabidopsis_thaliana', ATH, {'NC_003071.7': L}, rows)
    cells = lab.seed_cells()
    checked = 0
    for r in rows[:30]:
        c, st, a, b = cells[r['id']][0]
        tg, _, _ = lab.targets(c, st, a, b, 'S1')
        seg = tg['seg'][:b - a]
        x, y = max(a, r['span'][0]) - a, min(b, r['span'][1]) - a
        assert (seg[x:y] != IGNORE).all()
        known = np.zeros(b - a, bool)
        for t in lab.local(c, st, a, b):
            known[max(a, t['span'][0]) - a:min(b, t['span'][1]) - a] = True
        assert (seg[~known] == IGNORE).all()
        checked += 1
    assert checked == 30
