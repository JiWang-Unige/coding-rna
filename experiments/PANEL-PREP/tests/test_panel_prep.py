"""PANEL-PREP invariants (CPU, no data needed except the ont-only checks, which are skipped when the prepared files are absent)."""
import json

import pytest

import panel_common as P
import panel_exclusion as PX
import panel_exposure as PE
import panel_homology as PH


def rows(lengths, role='assembled-molecule', loc='Chromosome', unit='Primary Assembly', refseq=True):
    return [dict(name=str(i), role=role, location=loc, refseq=(f'NC_{i:06d}.1' if refseq else 'na'), genbank=f'CM{i:06d}.1', unit=unit, length=n) for i, n in enumerate(lengths)]


def test_main_sequence_rule_chromosome_branch():
    rep = rows([50_000_000, 30_000_000, 4_000_000]) + rows([60_000_000], loc='Chloroplast', unit='non-nuclear') + rows([12_000_000], role='unplaced-scaffold', loc='na')
    main, rule = P.main_sequences(rep)
    assert sorted(r['length'] for r in main) == [30_000_000, 50_000_000] and rule.startswith('assembled-molecule')     # 4 Mb, organelle and 12 Mb unplaced dropped


def test_main_sequence_rule_scaffold_fallback():
    rep = rows([88_000_000, 44_000_000, 9_999_999, 500_000], role='unplaced-scaffold', loc='na') + rows([400_000], role='unlocalized-scaffold', loc='Mitochondrion', unit='non-nuclear')
    main, rule = P.main_sequences(rep)
    assert sorted(r['length'] for r in main) == [44_000_000, 88_000_000] and rule.startswith('no assembled-molecule')  # only >= 10 Mb


def test_seq_name_refseq_vs_genbank():
    assert P.seq_name(rows([1])[0]) == 'NC_000000.1'
    assert P.seq_name(rows([1], refseq=False)[0]) == 'CM000000.1'


def test_species_registry_matches_m34_exclusion_list():
    sp = P.species()
    assert len(sp) == 12 and len({r['accession'] for r in sp.values()}) == 12
    px = {p[2]: p for p in PX.PANEL}
    for r in sp.values():
        assert r['accession'] in px and px[r['accession']][0] == r['genus'] and px[r['accession']][4] == 'main'
        assert r['ftp_base'].endswith(f"{r['accession']}_{r['assembly_name']}") and r['tier'] in ('TH', 'SH', 'AH')
    assert sum(r['tier'] == 'TH' for r in sp.values()) == 6 and sum(r['tier'] == 'SH' for r in sp.values()) == 5 and sum(r['tier'] == 'AH' for r in sp.values()) == 1


def test_exposure_audit_no_genus_or_accession_overlap_with_m34_training():
    a = PE.audit()
    assert not a['problems'] and not a['m34_direction']['violations']
    assert all(not v['m34_training_overlap']['genus'] for v in a['panel'].values())
    assert {k for k, v in a['panel'].items() if v['m34_training_overlap']['family']} == {'zizania_latifolia', 'themeda_triandra'}          # Poaceae only
    assert {k for k, v in a['panel'].items() if v['m30_adaptation_exposure']} == {'dryas_octopetala', 'eschscholzia_californica'}
    assert {k for k, v in a['panel'].items() if v['sealed_overlap']} == {'zizania_latifolia', 'themeda_triandra'}                         # Setaria: Poaceae family only
    assert set(a['panel']) == set(json.load(open(P.HERE / 'tool_exposure.json'))) - {'_legend'}


def test_homology_binning():
    assert [PH.binned(x) for x in (None, 0, 29.9, 30, 45, 69.9, 70, 89.9, 90, 100)] == ['no_hit', '<30', '<30', '30-40', '40-50', '50-70', '70-90', '70-90', '>=90', '>=90']
    s = PH.summarize(4, [95.0, 55.0, 31.0])
    assert s['frac']['no_hit'] == 0.25 and s['frac_ge']['50'] == 0.5 and s['frac_ge']['70'] == 0.25 and s['median_best_identity_with_hit'] == 55.0


@pytest.mark.skipif(not (P.PREP / 'dryas_octopetala' / 'splits.json').exists(), reason='prepared panel data not present')
def test_prepared_panel_consistent():
    for key, row in P.species().items():
        p = P.PREP / key / 'splits.json'
        if not p.exists():
            continue
        s = json.load(open(p))
        lens = json.load(open(P.PREP / key / 'lengths.json'))
        assert s['train'] == [] and s['val'] == [] and set(s['test']) == set(lens) == set(s['main_sequences'])
        assert all(lens[c] == s['main_sequences'][c]['length'] for c in lens) and s['test_bp'] == sum(lens.values())
        assert s['accession'] == row['accession'] and s['md5'][f"{s['assembly']}_genomic.gff.gz"]['status'] == 'OK'
        seen = set()
        for line in open(P.PREP / key / 'reference.gff3'):
            if line[0] != '#':
                f = line.split('\t')
                seen.add(f[0])
                assert int(f[4]) <= lens[f[0]]
        assert seen <= set(lens)
        # the FASTA headers are exactly the kept sequences
        heads = [l[1:].strip() for l in open(P.PREP / key / 'genome.fa') if l[0] == '>']
        assert heads == sorted(lens)
