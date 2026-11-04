"""M39 CPU checks of protocol 29 §5.3 that need only the prepared data and the label catalogs (no model):
  V1  prepared chromosomes vs protocol §3.2 (S / D / validation chromosome), raw md5 status;
  V2  catalog counts: primary coding genes, CDS-assessable references, eligible genes and eligible / assessable, non-standard splice fraction, CDS_partial;
      gate: eligible / assessable >= 90 % and non-standard splice <= 5 %;
  V3  Drosophila: coding genes parsed vs NCBI 13,962, sequence-structure failure reasons (read-through / non-ATG start etc. are not eligible);
  V4  NT-v2 window coverage of eligible genes: CDS span <= 12,288 / <= 24,064 / <= 49,152 / > 49,152 bp, genes with > 128 exons, introns < 20 bp;
      gate E6: span > 24,064 bp <= 15 %;
  V6  validation crops: CDS-assessable primary references of the catalog and of the frozen reference builder (as scored); gate E8: >= 400;
  V8  schedules: arms' labeled-gene counts, U and F share one schedule, no forbidden contig, species counts within 5 sigma of equal weights.
Writes <OUT>/checks/checks.json and prints one line per species.  usage: m39_checks.py
"""
import hashlib
import json
from collections import Counter

import numpy as np

import m39_common as C

PROTOCOL_HELDOUT = {      # protocol 29 §3.2, precomputed from the NCBI sequence reports on 2026-11-04
    'drosophila_melanogaster': ('NT_033778.4', 'NC_004354.4'), 'tribolium_castaneum': ('NC_087399.1', 'NC_087400.1'),
    'anopheles_gambiae': ('NC_064602.1', None), 'bombyx_mori': ('NC_085114.1', None), 'plodia_interpunctella': ('NC_071307.1', None),
    'spodoptera_frugiperda': ('NC_064215.1', None),                      # alternate (replaces bombyx_mori); held-out computed by the same rule on 11-04
    'apis_cerana': ('NC_083861.1', None), 'acromyrmex_echinatior': ('NC_141861.1', None), 'ips_typographus': ('NC_138430.1', None),
    'aethina_tumida': ('NC_065439.1', None), 'cloeon_dipterum': ('NC_088788.1', None)}


def introns(t):
    cds = sorted(t['CDS'])
    return [b[0] - a[1] for a, b in zip(cds, cds[1:])]


def species_checks(key, row, tab):
    sp = tab[key]
    splits = json.load(open(C.PREP / key / 'splits.json'))
    summ = json.load(open(C.CATALOG / key / 'summary.json'))
    cat = [json.loads(l) for l in open(C.CATALOG / key / 'transcripts.jsonl')]
    md5 = open(C.raw_dir(row) / 'local_md5.txt').read()
    held, d = PROTOCOL_HELDOUT[key]
    v1 = dict(heldout=splits['heldout_chromosome'], heldout_matches_protocol=splits['heldout_chromosome'] == held,
              d=splits['heldout_d'], d_matches_protocol=(splits['heldout_d'] == ([d] if d else [])), raw_md5_ok=md5.count(' OK') == 3,
              main_chromosomes=len(splits['main_chromosomes']), main_bp=splits['main_bp'], assembly_bp=splits['assembly_bp'],
              main_coverage=splits['main_bp'] / splits['assembly_bp'], train_contigs=len(splits['train']), train_bp=splits['train_bp'],
              validation_contig=splits['validation_contig'], validation_crop=splits['validation_crop'], discarded_bp=splits['discarded_bp'],
              dropped_trans_spliced_transcripts=len(splits.get('dropped_trans_spliced_transcripts', [])))
    ass = [t for t in cat if t['assessable']]
    el = [t for t in cat if t['eligible']]
    nonstd = sum('nonstandard_splice' in t['sequence_structure_failures'] for t in ass)
    v2 = dict(primary_coding=len(cat), assessable=len(ass), eligible=len(el), eligible_over_assessable=len(el) / len(ass),
              nonstandard_splice_fraction=nonstd / len(ass), cds_partial=sum(t['CDS_partial_attribute'] for t in cat),
              overlapping_primary=sum(t['overlapping_primary'] for t in cat),
              failure_reasons=dict(Counter(f for t in cat if t['assessable'] and not t['eligible'] for f in (t['sequence_structure_failures'] or ['overlap_or_partial']))),
              train_eligible=summ['stats']['train']['eligible_positive'], nested_labeled=summ['nested_labeled_counts'],
              n_coding_genes_parsed=summ['n_coding_genes_parsed'], ncbi_protein_coding=row['protein_coding'])
    v2['gate_ok'] = v2['eligible_over_assessable'] >= 0.90 and v2['nonstandard_splice_fraction'] <= 0.05
    spans = np.array([t['span'][1] - t['span'][0] for t in el])
    ints = [x for t in el for x in introns(t)]
    v4 = dict(eligible=len(el), span_median=float(np.median(spans)), span_p90=float(np.percentile(spans, 90)), span_max=int(spans.max()),
              le_12288=float((spans <= 12288).mean()), le_24064=float((spans <= 24064).mean()), le_49152=float((spans <= 49152).mean()),
              gt_24064=float((spans > 24064).mean()), gt_49152=float((spans > 49152).mean()),
              exons_gt_128=sum(t['n_exons'] > 128 for t in cat), max_exons_eligible=max(t['n_exons'] for t in el),
              introns=len(ints), introns_lt_20=sum(x < 20 for x in ints), intron_min=min(ints) if ints else None, intron_median=float(np.median(ints)) if ints else None)
    v4['gate_ok'] = v4['gt_24064'] <= 0.15
    out = dict(species=key, role=row['role'], V1=v1, V2=v2, V4=v4)
    if row['role'] == 'new':
        out['V6'] = dict(validation_contig=sp['val'][0], crop_bp=sp['lengths'][sp['val'][0]], catalog_assessable=summ['stats']['val']['cds_assessable_references'],
                         catalog_eligible=summ['stats']['val']['eligible_positive'])
    if key == 'drosophila_melanogaster':
        out['V3'] = dict(n_coding_genes_parsed=summ['n_coding_genes_parsed'], ncbi_protein_coding=row['protein_coding'],
                         structure_failures=dict(Counter(f for t in cat for f in t['sequence_structure_failures'])))
    return out


def reference_counts(tab):
    """Frozen reference builder on every evaluation contig (exactly what m39_val scores against)."""
    import m39_train_common as T
    MC = T.register()
    res = {}
    for key, sp in tab.items():
        for c in sp['val'] + sp['heldout_d']:
            ref = MC.build_reference(c)
            res[c] = dict(species=key, references=len(ref['refs']), length=MC.LENGTH[c], read_ok=len(MC.read_sequence(c)) == MC.LENGTH[c])
    return res


def schedule_checks():
    import m39_schedule as S
    out = {}
    sdir = C.OUT / 'schedules'
    summ = json.load(open(sdir / 'schedule_summary.json'))
    md5 = {a: hashlib.md5(open(sdir / f'schedule_{a}.jsonl', 'rb').read()).hexdigest() for a in S.ARMS}
    forbidden = C.forbidden_contigs()
    for arm in S.ARMS:
        rows = [json.loads(l) for l in open(sdir / f'schedule_{arm}.jsonl')]
        n, k = len(rows), len(S.ARMS[arm]['species'])
        cnt = Counter(r['species'] for r in rows)
        sd = (n * (1 / k) * (1 - 1 / k)) ** .5
        out[arm] = dict(md5=md5[arm], draws=n, unique_windows=summ[arm]['unique_windows'], labeled_genes_total=sum(summ[arm]['labeled_genes'].values()),
                        labeled_genes=summ[arm]['labeled_genes'], eligible_genes=summ[arm]['eligible_genes'],
                        windows_with_positive_fraction=summ[arm]['windows_with_positive_fraction'],
                        labeled_genes_exposed_unique=summ[arm]['labeled_genes_exposed_unique'],
                        forbidden_hits=sum(r['seqid'] in forbidden for r in rows), max_width=max(r['end'] - r['start'] for r in rows),
                        species_within_5sd=all(abs(c - n / k) <= 5 * sd for c in cnt.values()))
    out['U_F_identical_schedule'] = md5['I10_k50U'] == md5['I10_k50F']
    return out


def main():
    tab = C.species_table()
    rows = C.roster()
    res = dict(species={k: species_checks(k, rows[k], tab) for k in rows})
    res['references'] = reference_counts(tab)
    for k, r in res['species'].items():
        if 'V6' in r:
            r['V6']['scored_references'] = res['references'][r['V6']['validation_contig']]['references']
            r['V6']['gate_ok'] = r['V6']['scored_references'] >= C.VAL_MIN_REFS
    res['schedules'] = schedule_checks()
    (C.OUT / 'checks').mkdir(parents=True, exist_ok=True)
    json.dump(res, open(C.OUT / 'checks' / 'checks.json', 'w'), indent=1)
    for k, r in res['species'].items():
        v2, v4 = r['V2'], r['V4']
        print(json.dumps(dict(species=k, heldout_ok=r['V1']['heldout_matches_protocol'] and r['V1']['d_matches_protocol'], md5_ok=r['V1']['raw_md5_ok'],
                              assessable=v2['assessable'], eligible=v2['eligible'], elig_frac=round(v2['eligible_over_assessable'], 4),
                              nonstd=round(v2['nonstandard_splice_fraction'], 4), V2_ok=v2['gate_ok'], gt24064=round(v4['gt_24064'], 4), V4_ok=v4['gate_ok'],
                              exons_gt_128=v4['exons_gt_128'], introns_lt_20=v4['introns_lt_20'],
                              val_refs=r.get('V6', {}).get('scored_references'), V6_ok=r.get('V6', {}).get('gate_ok'))))
    print(json.dumps(dict(schedules={a: {k: v for k, v in s.items() if k in ('labeled_genes_total', 'unique_windows', 'forbidden_hits', 'species_within_5sd', 'md5')}
                                     for a, s in res['schedules'].items() if isinstance(s, dict)}, U_F_identical=res['schedules']['U_F_identical_schedule'])))


if __name__ == '__main__':
    main()
