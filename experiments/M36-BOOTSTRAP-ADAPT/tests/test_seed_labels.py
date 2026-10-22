"""Seed-GFF -> M34 label interface (m36_seed_labels).  Synthetic part (no data): GFF round trip, unknown mask, seed roles vs GFF, F1 has no IGNORE,
pseudo-background overlay, S3 merge.  Data part (skipped without the M34 prepared Brachypodium data): F1 targets of a SeedLabeler fed with the
catalog genes of a window equal `m34_labels.Labeler.targets(..., 'FULL')` exactly; S-arm targets equal FULL on the seed loci and are IGNORE elsewhere.
Run: pytest -q test_seed_labels.py"""
import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

import m36_common as W
import m36_seed_labels as SL
import manifest as MF
from src.m28.labels import IGNORE, oriented_parts

SAFE = ['GCT', 'GAA', 'AAA', 'CCT', 'GGC', 'TTC', 'GAC', 'CAA']


def make_gene(rng, n_exons):
    """-> (oriented gene DNA string, exon pieces as (start, end) relative to the gene start, phases) : legal literal-DNA complete CDS"""
    n_codons = int(rng.integers(12, 30))
    cds = 'ATG' + ''.join(rng.choice(SAFE, n_codons)) + 'TAA'
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


def synthetic_contig(rng, length=60000, n_genes=6):
    """random contig with planted genes on both strands -> (dna, seed rows); the '-' genes are reverse-complemented into place"""
    dna = list(rng.choice(list('ACGT'), length))
    rows, cursor = [], 2000
    for k in range(n_genes):
        g, exons = make_gene(rng, int(rng.integers(1, 6)))
        strand = '+' if k % 2 == 0 else '-'
        start = cursor
        if strand == '+':
            dna[start:start + len(g)] = list(g)
            cds = [(start + a, start + b) for a, b in exons]
        else:
            rc = MF.rc(g)
            dna[start:start + len(g)] = list(rc)
            cds = sorted((start + len(g) - b, start + len(g) - a) for a, b in exons)
        rows.append(SL.seed_row('SYN', strand, SL.phases(cds, strand), f'g{k}', identity=0.9, coverage=1.0, source_protein=f'P{k}', source_species='x'))
        cursor = start + len(g) + int(rng.integers(500, 3000))
    dna = ''.join(dna)
    for r in rows:
        failures, _ = MF.chain_checks(r, dna)
        assert not failures, (r['id'], failures)
    return dna, rows


@pytest.fixture(scope='module')
def syn():
    rng = np.random.default_rng(20261021)
    dna, rows = synthetic_contig(rng)
    td = tempfile.mkdtemp()
    fa = Path(td) / 'genome.fa'
    with open(fa, 'w') as f:
        f.write('>SYN\n')
        for i in range(0, len(dna), 80):
            f.write(dna[i:i + 80] + '\n')
    return dict(dna=dna, rows=rows, fa=fa, lengths={'SYN': len(dna)}, td=Path(td))


def write_gff(rows, path):
    with open(path, 'w') as f:
        f.write('##gff-version 3\n')
        for r in rows:
            a, b = r['span']
            common = f"source_protein={r['source_protein']};source_species={r['source_species']};identity={r['identity']};coverage={r['coverage']}"
            f.write(f"SYN\tx\tgene\t{a + 1}\t{b}\t.\t{r['strand']}\t.\tID={r['id']};{common}\n")
            f.write(f"SYN\tx\tmRNA\t{a + 1}\t{b}\t.\t{r['strand']}\t.\tID={r['id']}.m;Parent={r['id']};{common}\n")
            for j, (x, y, p) in enumerate(r['CDS']):
                f.write(f"SYN\tx\tCDS\t{x + 1}\t{y}\t.\t{r['strand']}\t{p}\tParent={r['id']}.m\n")


def test_gff_round_trip(syn):
    p = syn['td'] / 'seeds.gff3'
    write_gff(syn['rows'], p)
    back = SL.load_seed_gff(p)
    assert [(r['id'], r['seqid'], r['strand'], r['CDS']) for r in back] == [(r['id'], r['seqid'], r['strand'], r['CDS']) for r in syn['rows']]
    json.dump(syn['rows'], open(syn['td'] / 'seeds.json', 'w'))
    back2 = SL.load_seed_json(syn['td'] / 'seeds.json')
    assert [r['CDS'] for r in back2] == [r['CDS'] for r in syn['rows']]


def expected_roles(row, win):
    """roles implied by the GFF chain (the frozen label semantics): CDS phases, donor/acceptor on the last/first exon base, start/stop, intron phases"""
    L, n = win['width'], win['valid_bases']
    parts = oriented_parts(row, win)
    seg = np.zeros(L, dtype=np.int64); region = np.zeros(L, dtype=np.int64)
    off = 0
    for j, (x, y) in enumerate(parts):
        seg[x:y] = 1 + (off + np.arange(y - x)) % 3; region[x:y] = 1
        off += y - x
        if j < len(parts) - 1:
            z = parts[j + 1][0]
            seg[y:z] = 4 + (off - 1) % 3; region[y:z] = 2
    off = 0
    for (x, y), (z, _) in zip(parts, parts[1:]):
        off += y - x
        seg[y - 1] = 7 + (off - 1) % 3; seg[z] = 10 + off % 3
    a, b = parts[0][0], parts[-1][1]
    seg[a:a + 3] = 13; seg[b - 3:b] = 14
    return seg, region, (a, b)


def test_seed_roles_and_unknown_mask(syn):
    lab = SL.SeedLabeler('syn', syn['fa'], syn['lengths'], syn['rows'])
    n_checked = 0
    for strand in '+-':
        for a, b in lab.grid('SYN'):
            local = lab.local('SYN', strand, a, b)
            if not local:
                continue
            (tg, win, core) = lab.targets('SYN', strand, a, b, 'S2')
            n = b - a
            known = np.zeros(win['width'], dtype=bool)
            for t in local:
                seg_e, reg_e, (x, y) = expected_roles(t, win)
                x, y = max(0, x), min(n, y)
                if y <= x:
                    continue
                assert np.array_equal(tg['seg'][x:y], seg_e[x:y]), t['id']          # seed labels == GFF chain
                assert np.array_equal(tg['region'][x:y], reg_e[x:y]), t['id']
                known[x:y] = True
                n_checked += 1
            assert (tg['seg'][~known] == IGNORE).all() and (tg['region'][~known] == IGNORE).all()      # everything else unknown
            assert (tg['seg'][n:] == IGNORE).all()
            assert (tg['seg'][known] >= 0).all()                                                       # seed bases are never unknown (ACGT synthetic)
    assert n_checked >= len(syn['rows'])


def test_F1_has_no_ignore_and_matches_S_on_seeds(syn):
    lab = SL.SeedLabeler('syn', syn['fa'], syn['lengths'], syn['rows'])
    for strand in '+-':
        for a, b in lab.grid('SYN'):
            n = b - a
            F, _, _ = lab.targets('SYN', strand, a, b, 'F1')
            S, _, _ = lab.targets('SYN', strand, a, b, 'S1')
            assert (F['seg'][:n] >= 0).all() and (F['region'][:n] >= 0).all()           # F1: no unknown anywhere
            assert (F['seg'][n:] == IGNORE).all()
            m = S['seg'] >= 0
            assert np.array_equal(F['seg'][m], S['seg'][m]) and np.array_equal(F['region'][m], S['region'][m])
            assert (F['seg'][~m][:max(0, n - m.sum())] == 0).all() if n > m.sum() else True   # F1 fills the unknown positions with background


def test_pseudo_background_overlay(syn):
    rows = syn['rows']
    g = rows[0]                                   # '+' seed
    L = syn['lengths']['SYN']
    bg = {('SYN', '+'): [(0, 1000), (g['span'][1] + 100, g['span'][1] + 600)], ('SYN', '-'): [(5000, 5200)]}
    lab = SL.SeedLabeler('syn', syn['fa'], syn['lengths'], rows, background=bg)
    a, b = lab.grid('SYN')[0]
    S, win, _ = lab.targets('SYN', '+', a, b, 'S1')
    assert (S['seg'][0:1000] == 0).all() and (S['region'][0:1000] == 0).all()
    x, y = g['span'][1] + 100 - a, g['span'][1] + 600 - a
    assert (S['seg'][x:y] == 0).all()
    assert (S['seg'][1000:1200] == IGNORE).all()
    Sm, _, _ = lab.targets('SYN', '-', a, b, 'S1')
    x, y = MF.orient_interval(5000, 5200, a, b, '-')
    assert (Sm['seg'][x:y] == 0).all() and (Sm['seg'][y:y + 50] == IGNORE).all()
    # pseudo-background from a probability vector on the '-' strand lands at the right genomic place
    p = np.zeros(L); p[L - 5200:L - 5000] = 1.0
    assert SL.pseudo_background(p, '-', L) == [(5000, 5200)]


def test_merge_S3(syn):
    rows = syn['rows']
    s2 = rows[:3]
    s1 = [dict(rows[1], id='alt1'), dict(rows[4], id='alt4')]      # alt1 conflicts with an S2 seed, alt4 does not
    merged, dropped = SL.merge_seed_sets(s2, s1)
    assert dropped == 1 and {r['id'] for r in merged} == {'g0', 'g1', 'g2', 'alt4'}
    assert SL.arm_seeds('S3', s1_rows=s1, s2_rows=s2)[1] == 1
    assert [r['id'] for r in SL.arm_seeds('F1', s1_rows=s1, s2_rows=s2)[0]] == ['alt1', 'alt4']


# ---- data part: equality with the M34 FULL labels -------------------------------------------------------------------------------------
import m34_common as C
HAVE = (C.PREP / 'brachypodium_distachyon' / 'splits.json').exists()


@pytest.mark.skipif(not HAVE, reason='prepared M34 data not present')
def test_equals_m34_full_on_brachypodium():
    import m34_labels as LB
    key = 'brachypodium_distachyon'
    sp = C.species_table()[key]
    lab34 = LB.Labeler(key, sp['val'], 1.0)
    contig = sp['val'][0]
    seeds = [SL.seed_row(t['seqid'], t['strand'], t['CDS'], t['id'], gene_id=t['gene_id']) for t in lab34.catalog.values() if t['eligible']]
    lab36 = SL.SeedLabeler(key, sp['dir'] / 'genome.fa', {contig: sp['lengths'][contig]}, seeds)
    checked = {'+': 0, '-': 0}
    for strand in '+-':
        for a, b in lab34.grid(contig):
            local = lab34.local(contig, strand, a, b)
            if not local or not all(t['eligible'] for t in local) or checked[strand] >= 4:
                continue
            FULL, win34, core34 = lab34.targets(contig, strand, a, b, 'FULL')
            F1, win36, core36 = lab36.targets(contig, strand, a, b, 'F1')
            assert core34 == core36
            assert np.array_equal(FULL['seg'], F1['seg']) and np.array_equal(FULL['region'], F1['region']), (contig, strand, a)
            S, _, _ = lab36.targets(contig, strand, a, b, 'S2')
            known = np.zeros(win36['width'], dtype=bool)
            for t in local:
                x, y = MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, strand)
                known[x:y] = True
            assert np.array_equal(S['seg'][known], FULL['seg'][known]) and (S['seg'][~known] == IGNORE).all() and (S['region'][~known] == IGNORE).all()
            checked[strand] += 1
    assert checked['+'] >= 1 and checked['-'] >= 1, checked
    print('windows compared', checked)
