"""Exposure audit of the 12 panel species against OUR training-side sets (the mirror image of M34's panel_exclusion audit):
genus / family / accession overlap with the M34 10 training species, the M30 adaptation species (Dryas, Eschscholzia) and the project-sealed Setaria;
plus consistency with M34's exclusion list (panel_exclusion.PANEL) and the per-tool exposure codes transcribed from doc 10 §1.7 (tool_exposure.json).
usage: panel_exposure.py [out.json]      (pure Python, no data needed)
"""
import json
import sys

import panel_common as P
import panel_exclusion as PX


def audit():
    ex = P.exposure_sets()
    side = [(grp, e) for grp, lst in ex.items() for e in lst]
    px = {p[2]: p for p in PX.PANEL if p[2]}
    tools = json.load(open(P.HERE / 'tool_exposure.json'))
    out, problems = {}, []
    for key, row in P.species().items():
        hits = dict(same_accession=[], same_genus=[], same_family=[])
        for grp, e in side:
            tag = f"{e['species']} [{e['role']}]"
            if e['accession'] and e['accession'] == row['accession']:
                hits['same_accession'].append(tag)
            if e['genus'] == row['genus']:
                hits['same_genus'].append(tag)
            if e['family'] == row['family']:
                hits['same_family'].append(tag)
        m34_hits = [t for t in hits['same_genus'] + hits['same_accession'] if 'M34 training' in t]
        if m34_hits:
            problems.append((key, m34_hits))
        in_px = px.get(row['accession'])
        if not in_px or in_px[0] != row['genus'] or in_px[4] != 'main':
            problems.append((key, 'not in panel_exclusion.PANEL as a main-panel entry with the same genus'))
        out[key] = dict(organism=row['organism'], accession=row['accession'], genus=row['genus'], family=row['family'], tier=row['tier'], tier_note=row['tier_note'],
                        m34_training_overlap=dict(genus=bool(m34_hits), family=[t for t in hits['same_family'] if 'M34 training' in t]),
                        m30_adaptation_exposure=[t for t in hits['same_accession'] if 'M30' in t],
                        sealed_overlap=[t for t in hits['same_genus'] + hits['same_family'] if 'sealed' in t],
                        hits=hits, in_m34_exclusion_list=bool(in_px), tool_exposure=tools[key])
    # the M34 direction (training -> panel) must still hold
    training = [(e['genus'], e['accession'], e['family']) for e in ex['m34_training']]
    m34 = PX.audit(training)
    return dict(panel=out, m34_direction=dict(violations=m34['violations'], family_overlaps=sorted({(x['training_genus'], x['panel_genus']) for x in m34['family_overlaps'] if x['panel_tier'] == 'main'})),
                exposure_sets=ex, problems=problems,
                rule='no panel species may share a genus or an accession with an M34 training species (family overlaps are reported, not forbidden); '
                     'M30 adaptation exposure (Dryas, Eschscholzia) is flagged: any checkpoint that saw their label-set-A chromosomes must not be scored on those chromosomes')


if __name__ == '__main__':
    a = audit()
    out = sys.argv[1] if len(sys.argv) > 1 else str(P.HERE.parents[0] / 'results' / 'exposure_audit.json')
    json.dump(a, open(out, 'w'), indent=1)
    for k, v in a['panel'].items():
        print(f"{k:26s} {v['tier']:3s} genus_overlap={v['m34_training_overlap']['genus']} family={v['m34_training_overlap']['family']} m30={v['m30_adaptation_exposure']}")
    print('problems:', a['problems'])
    assert not a['problems'], a['problems']
