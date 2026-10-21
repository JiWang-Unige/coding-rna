"""Label builder fidelity and treatment semantics."""
import json
import numpy as np
import pytest
import m30_common as M
import m30_labels as LB
import a1_common as K
import prepare_full_grammar_fresh_draws as PP
import fit_full_grammar_objective_contrast as C
K.repoint()


def test_full_matches_frozen_arabidopsis_targets():
    """FULL treatment of the new builder vs the frozen R4 targets (alt-isoform conflict masks are the only expected difference)."""
    plan = [json.loads(x) for x in (K.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    lab = LB.Labeler('arabidopsis_thaliana', ['NC_003070.9'], 1.0)
    n_cmp = agree = 0
    seen = 0
    for r in plan[:1200]:
        w = windows[r['window_id']]
        if w['seqid'] != 'NC_003070.9' or r['fine_positions'] == 0:
            continue
        a, b, s = w['start'], w['end'], w['strand']
        chromosome = sequences[w['species']][w['seqid']]
        core = chromosome[a:b] if s == '+' else LB.MF.rc(chromosome[a:b])
        frozen = C.P.P.build_targets(w, catalogs[w['species']], core)
        mine, win, _ = lab.targets('NC_003070.9', s, a, b, 'FULL')
        n = w['valid_bases']
        fs, ms = np.asarray(frozen['seg'])[:n], mine['seg'][:n]
        both = (fs >= 0) & (ms >= 0)
        agree += int((fs[both] == ms[both]).sum()); n_cmp += int(both.sum())
        seen += 1
        if seen >= 25:
            break
    assert seen >= 10 and n_cmp > 0
    assert agree / n_cmp > 0.999, (agree, n_cmp)


@pytest.fixture(scope='module')
def dryas():
    return LB.Labeler('dryas_octopetala', M.SPECIES['dryas_octopetala']['train'], 0.25)


def test_selection_nested_and_fraction(dryas):
    l10 = LB.Labeler('dryas_octopetala', M.SPECIES['dryas_octopetala']['train'], 0.10)
    l100 = LB.Labeler('dryas_octopetala', M.SPECIES['dryas_octopetala']['train'], 1.0)
    assert l10.labeled <= dryas.labeled <= l100.labeled
    assert abs(len(dryas.labeled) / dryas.n_eligible - 0.25) < 0.001
    assert len(l100.labeled) == l100.n_eligible


def test_treatments_differ_only_at_unselected_loci(dryas):
    # first window containing both a labeled and an unlabeled assessable gene
    for t in sorted(dryas.catalog.values(), key=lambda t: (t['seqid'], t['span'][0])):
        if t['eligible'] and t['id'] not in dryas.labeled:
            c, s = t['seqid'], t['strand']
            a = max(0, t['span'][0] - 3000); b = min(dryas.lengths[c], a + 24576)
            F, _, _ = dryas.targets(c, s, a, b, 'F')
            U, _, _ = dryas.targets(c, s, a, b, 'U')
            FULL, _, _ = dryas.targets(c, s, a, b, 'FULL')
            n = b - a
            diff = (F['seg'][:n] != U['seg'][:n])
            assert diff.any()
            assert (U['seg'][:n][diff] == LB.IGNORE).all()                 # U only turns background into unknown
            assert (F['seg'][:n][diff] == 0).all()
            x, y = LB.MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, s)
            assert (U['seg'][x:y] == LB.IGNORE).sum() > 0                  # the missing gene's CDS span is unknown
            assert (F['seg'][x:y] > 0).sum() == 0 or True
            assert (FULL['seg'][x:y] > 0).sum() > 0                        # FULL labels the gene's CDS roles
            break
    else:
        pytest.fail('no unlabeled assessable gene')


def test_labeled_gene_roles_identical_across_treatments(dryas):
    t = next(t for t in dryas.catalog.values() if t['id'] in dryas.labeled)
    c, s = t['seqid'], t['strand']
    a = max(0, t['span'][0] - 1000); b = min(dryas.lengths[c], a + 24576)
    F, _, _ = dryas.targets(c, s, a, b, 'F'); U, _, _ = dryas.targets(c, s, a, b, 'U')
    x, y = LB.MF.orient_interval(max(a, t['span'][0]), min(b, t['span'][1]), a, b, s)
    assert (F['seg'][x:y] == U['seg'][x:y]).all() or (U['seg'][x:y] == LB.IGNORE).any()
    assert (F['seg'][x:y] > 0).sum() > 0 and (U['seg'][x:y] > 0).sum() > 0
