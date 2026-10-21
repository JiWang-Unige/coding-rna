"""Independent-panel exclusion list (docs/review_20261016/10_independent_plant_panel.md §2.2/2.3) + sealed/other project species.
No M34 training species may share an accession or a GENUS with any entry (family overlaps are reported, not forbidden)."""
PANEL = [  # (genus, species_epithet_or_None, accession, family, tier)
    ('Dryas', 'octopetala', 'GCF_963921425.1', 'Rosaceae', 'main'), ('Eschscholzia', 'californica', 'GCF_041429995.1', 'Papaveraceae', 'main'),
    ('Akebia', 'trifoliata', 'GCF_017979445.1', 'Lardizabalaceae', 'main'), ('Hakea', 'prostrata', 'GCF_041222555.1', 'Proteaceae', 'main'),
    ('Chenopodium', 'album', 'GCF_948465745.1', 'Amaranthaceae', 'main'), ('Zizania', 'latifolia', 'GCA_048537475.1', 'Poaceae', 'main'),
    ('Themeda', 'triandra', 'GCA_057662145.1', 'Poaceae', 'main'), ('Juncus', 'effusus', 'GCA_027726005.1', 'Juncaceae', 'main'),
    ('Rhynchospora', 'breviuscula', 'GCA_027562975.1', 'Cyperaceae', 'main'), ('Acorus', 'gramineus', 'GCA_030737835.1', 'Acoraceae', 'main'),
    ('Dipteris', 'shenzhenensis', 'GCF_047967505.1', 'Dipteridaceae', 'main'), ('Amborella', 'trichopoda', 'GCF_054954075.1', 'Amborellaceae', 'main'),
    ('Coptis', 'chinensis', 'GCF_015680905.1', 'Ranunculaceae', 'alt'), ('Nelumbo', 'lutea', 'GCF_056320005.1', 'Nelumbonaceae', 'alt'),
    ('Medicago', 'sativa', 'GCF_051527215.1', 'Fabaceae', 'alt'), ('Azolla', 'caroliniana', 'GCF_043228905.1', 'Salviniaceae', 'alt'),
    ('Adiantum', 'capillus-veneris', 'GCF_014529385.1', 'Pteridaceae', 'alt'), ('Lygodium', 'microphyllum', 'GCF_053477565.2', 'Lygodiaceae', 'alt'),
    ('Cibotium', 'barometz', 'GCF_041138505.1', 'Cibotiaceae', 'alt'), ('Canna', 'indica', 'GCA_034359265.1', 'Cannaceae', 'alt'),
    ('Stemona', 'tuberosa', 'GCA_055697745.1', 'Stemonaceae', 'alt'), ('Paspalum', 'saurae', 'GCA_036689595.1', 'Poaceae', 'alt'),
    ('Aristolochia', 'fimbriata', 'GCA_019845555.1', 'Aristolochiaceae', 'alt'), ('Nymphaea', 'thermarum', 'GCA_011799765.1', 'Nymphaeaceae', 'alt'),
    ('Sphagnum', 'compactum', 'GCA_965153455.1', 'Sphagnaceae', 'alt'), ('Ceratodon', 'purpureus', 'GCA_014871385.1', 'Ditrichaceae', 'alt'),
    ('Edaphochlamys', 'debaryana', 'GCA_016858145.1', 'Chlamydomonadaceae', 'alt'),
    ('Lactuca', None, None, 'Asteraceae', 'excluded-too-recent'),            # 2026-02 RefSeq Lactuca (doc 10 §2.1)
    ('Setaria', None, None, 'Poaceae', 'project-sealed'),                      # sealed project data: never used (M28-M33 rule)
]


def audit(training):
    """training: iterable of (genus, accession, family) -> dict with violations (must be empty) and family overlaps (reported)."""
    genera = {p[0] for p in PANEL}
    accs = {p[2] for p in PANEL if p[2]}
    viol, fam = [], []
    for g, a, f in training:
        if g in genera:
            viol.append(('genus', g))
        if a in accs:
            viol.append(('accession', a))
        for p in PANEL:
            if p[3] == f:
                fam.append(dict(training_genus=g, family=f, panel_genus=p[0], panel_tier=p[4]))
    return dict(violations=viol, family_overlaps=fam)
