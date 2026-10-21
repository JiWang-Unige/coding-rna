"""M34 data-preparation invariants (CPU).  Tests that need prepared data/catalogs on ont are skipped if those files are absent."""
import bisect
import json

import pytest

import m34_common as C
import panel_exclusion as PX


def rows(lengths, role='assembled-molecule', loc='Chromosome', unit='Primary Assembly'):
    return [dict(name=str(i), role=role, location=loc, refseq=f'NC_{i:06d}.1', unit=unit, length=n) for i, n in enumerate(lengths)]


def test_main_chromosome_rule_and_heldout_median():
    rep = rows([50_000_000, 30_000_000, 20_000_000, 4_000_000, 10_000_000]) + rows([60_000_000], loc='Chloroplast', unit='non-nuclear') + rows([9_000_000], role='unplaced-scaffold', loc='na')
    chroms = C.main_chromosomes(rep)
    assert sorted(r['length'] for r in chroms) == [10_000_000, 20_000_000, 30_000_000, 50_000_000]       # 4 Mb, organelle, unplaced dropped
    assert C.choose_heldout(chroms)['length'] == 20_000_000                                              # lower median of 4
    assert C.choose_heldout(chroms[:3])['length'] in (20_000_000, 30_000_000, 50_000_000)
    assert C.choose_heldout(rows([7_000_000]))['length'] == 7_000_000


def test_nested_labeled_selection():
    t = [dict(id=f'g{i}', seqid='c1' if i % 2 else 'c2', span=[i * 1000, i * 1000 + 500]) for i in range(1003)]
    s10, s25, s50, s100 = (C.select_labeled(t, k) for k in (.1, .25, .5, 1.0))
    assert s10 <= s25 <= s50 <= s100 and len(s100) == 1003
    assert len(s25) == -(-1003 // 4) and len(s50) == -(-1003 // 2)
    assert C.select_labeled(t, .25) == s25                                  # deterministic


def test_gene_cells_match_m30_rule():
    L = 200_000
    cells = C.grid(L)
    starts = [a for a, _ in cells]
    for span in ([5, 400], [12_000, 13_000], [24_500, 40_000], [150_000, 199_999], [0, 199_999], [100_000, 190_000]):
        i = bisect.bisect_right(starts, span[0]) - 1
        ref = [cells[j] for j in range(max(0, i - 2), min(len(cells), i + 3)) if cells[j][0] < span[1] and cells[j][1] > span[0]]
        assert C.gene_cells(span, L) == ref and ref


def test_panel_exclusion_audit():
    ok = PX.audit([('Populus', 'GCF_000002775.5', 'Salicaceae'), ('Brachypodium', 'GCF_000005505.3', 'Poaceae')])
    assert not ok['violations'] and any(x['panel_genus'] == 'Themeda' for x in ok['family_overlaps'])
    bad = PX.audit([('Dryas', 'GCF_000000000.1', 'Rosaceae'), ('Zea', 'GCA_057662145.1', 'Poaceae'), ('Setaria', 'x', 'Poaceae')])
    assert {v[1] for v in bad['violations']} == {'Dryas', 'GCA_057662145.1', 'Setaria'}


def test_species_json_is_clean():
    rows_ = C.new_species()
    assert len(rows_) == 8
    training = [(r['organism'].split()[0], r['accession'], C.FAMILY[k]) for k, r in rows_.items()]
    training += [(v['organism'].split()[0], v['accession'], v['family']) for v in C.LEGACY.values()]
    assert not PX.audit(training)['violations']
    assert len({r['accession'] for r in rows_.values()}) == 8


@pytest.mark.skipif(not (C.PREP / 'panel_exclusion_audit.json').exists(), reason='prepared data not present')
def test_prepared_splits_consistent():
    for key in C.new_species():
        p = C.PREP / key / 'splits.json'
        if not p.exists():
            continue
        s = json.load(open(p))
        assert s['heldout_chromosome'] not in s['train'] and s['validation_contig'] not in s['train']
        assert s['lengths'][s['validation_contig']] == min(C.VAL_BP, s['heldout_length'])
        assert set(s['lengths']) == set(s['train']) | {s['validation_contig']}
        seen = set()
        for line in open(C.PREP / key / 'reference.gff3'):
            if line[0] != '#':
                seen.add(line.split('\t', 1)[0])
        assert seen <= set(s['lengths']) and s['validation_contig'] in seen
        # no gene of the validation contig crosses the crop boundary
        for line in open(C.PREP / key / 'reference.gff3'):
            f = line.split('\t')
            if line[0] != '#' and f[0] == s['validation_contig'] and f[2] == 'gene':
                assert int(f[4]) <= s['lengths'][s['validation_contig']]


@pytest.mark.skipif(not (C.OUT / 'schedules' / 'schedule_summary.json').exists(), reason='schedules not generated')
def test_schedules_clean_and_deterministic():
    summ = json.load(open(C.OUT / 'schedules' / 'schedule_summary.json'))
    forbidden = set(C.DEV_CONTIGS) | {c for s in C.LEGACY.values() for c in s['val'] + s['never_train']}
    for arm, s in summ.items():
        assert s['draws'] == s['steps'] * s['batch']
        n = len(s['species'])
        for sp, d in s['draws_per_species'].items():
            assert abs(d / s['draws'] - 1 / n) < 0.02, (arm, sp)
        for line in open(C.OUT / 'schedules' / f'schedule_{arm}.jsonl'):
            r = json.loads(line)
            assert r['seqid'] not in forbidden and '__val' not in r['seqid']
