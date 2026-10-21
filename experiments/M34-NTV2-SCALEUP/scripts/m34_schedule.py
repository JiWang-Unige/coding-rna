"""M34 deterministic draw schedules for every arm (CPU only; reads data/m34_train/catalog/* written by m34_labels.py; no model).
Sampler law (one law for all arms, seed 1): per draw, species by the arm's weights (equal per species); then with probability 0.5 a window uniform over
that species' (training contig, strand, 24,576-bp grid window) pool, else 'gene-first': a uniform labeled eligible gene of that species, then a uniform grid
window among the (<=5) windows covering its CDS span, on the gene's strand (identical to the M30 sampler).  Validation/DEV/never-used contigs are not in the pool.
Outputs (ont ~/coding-rna/outputs/M34-NTV2-SCALEUP/schedules/): schedule_<arm>.jsonl, schedule_summary.json
usage: m34_schedule.py [arm ...]
"""
import bisect
import json
import sys
from collections import Counter

import numpy as np

import m34_common as C

STEPS, BATCH, SAMPLING_SEED, UNIFORM_FRACTION = 3072, 8, 1, 0.5
LEGACY = ['arabidopsis_thaliana', 'oryza_sativa']
NEW = ['populus_trichocarpa', 'brachypodium_distachyon', 'solanum_lycopersicum', 'phoenix_dactylifera', 'phaseolus_vulgaris', 'citrus_sinensis',
       'daucus_carota', 'dioscorea_rotundata']                                   # nested adding order, fixed in the protocol
ARMS = {
    'T2': dict(species=LEGACY, new_k=1.0, treatment='FULL'),
    'T4': dict(species=LEGACY + NEW[:2], new_k=1.0, treatment='FULL'),
    'T6': dict(species=LEGACY + NEW[:4], new_k=1.0, treatment='FULL'),
    'T10': dict(species=LEGACY + NEW, new_k=1.0, treatment='FULL'),
    'T10_k50U': dict(species=LEGACY + NEW, new_k=0.5, treatment='U'),
    'T10_k25U': dict(species=LEGACY + NEW, new_k=0.25, treatment='U'),
    'T10_k50F': dict(species=LEGACY + NEW, new_k=0.5, treatment='F'),
}


def load_species(key, k):
    """-> dict(contigs, lengths, genes=[(id, seqid, strand, span)] labeled at fraction k (training contigs only), cells per gene)."""
    sp = C.species_table()[key]
    tr = set(sp['train'])
    elig = [t for t in map(json.loads, open(C.CATALOG / key / 'transcripts.jsonl')) if t['eligible'] and t['seqid'] in tr]
    labeled = sorted(C.select_labeled(elig, k)) if k < 1 else sorted(t['id'] for t in elig)
    byid = {t['id']: t for t in elig}
    genes = [(i, byid[i]['seqid'], byid[i]['strand'], byid[i]['span']) for i in labeled]
    contigs = list(sp['train'])
    nwin = [len(C.grid(sp['lengths'][c])) for c in contigs]
    idx = {}                                                                 # (contig, strand) -> sorted [(start, end, id)] of labeled eligible genes
    for gid, c, st, span in genes:
        idx.setdefault((c, st), []).append((span[0], span[1], gid))
    for v in idx.values():
        v.sort()
    return dict(contigs=contigs, lengths=sp['lengths'], genes=genes, nwin=nwin, eligible=len(elig), index={k: ([x[0] for x in v], v) for k, v in idx.items()})


def positives(P, seqid, strand, a, b):
    """labeled eligible genes fully inside [a, b) with the 256-bp flanks (same rule as the manifest `positive_ids`)."""
    starts, v = P['index'].get((seqid, strand), ([], []))
    out = []
    for i in range(bisect.bisect_left(starts, a + C.FLANK), len(v)):
        s, e, gid = v[i]
        if s > b - C.FLANK:
            break
        if e <= b - C.FLANK:
            out.append(gid)
    return out


def draw_schedule(arm, spec, pools):
    rng = np.random.default_rng(SAMPLING_SEED)
    names = spec['species']
    w = np.full(len(names), 1.0 / len(names))
    rows = []
    total = STEPS * BATCH
    sp_idx = rng.choice(len(names), size=total, p=w)
    for i in range(total):
        key = names[int(sp_idx[i])]
        P = pools[key]
        mode = 'U'
        if rng.random() < UNIFORM_FRACTION or not P['genes']:
            tot = 2 * sum(P['nwin'])
            u = int(rng.integers(tot))
            strand = '+' if u < tot // 2 else '-'
            u %= tot // 2
            for c, n in zip(P['contigs'], P['nwin']):
                if u < n:
                    a, b = C.grid(P['lengths'][c])[u]
                    break
                u -= n
        else:
            mode = 'G'
            gid, c, strand, span = P['genes'][int(rng.integers(len(P['genes'])))]
            cells = C.gene_cells(span, P['lengths'][c])
            a, b = cells[int(rng.integers(len(cells)))]
        rows.append(dict(draw=i, step=i // BATCH + 1, species=key, seqid=c, strand=strand, start=a, end=b, mode=mode,
                         id=f'{key}|{c}|{a}|{strand}', **({'gene': gid} if mode == 'G' else {})))
    return rows


def main():
    out = C.OUT / 'schedules'
    out.mkdir(parents=True, exist_ok=True)
    arms = sys.argv[1:] or list(ARMS)
    sp_path = out / 'schedule_summary.json'
    summary = json.load(open(sp_path)) if sp_path.exists() else {}      # merge: arms can be generated in separate invocations
    cache = {}
    forbidden = set(C.DEV_CONTIGS) | {c for s in C.LEGACY.values() for c in s['val'] + s['never_train']} | \
        {s['val'][0] for s in C.species_table().values() if s['val']}
    for arm in arms:
        spec = ARMS[arm]
        pools = {}
        for key in spec['species']:
            k = 1.0 if key in LEGACY else spec['new_k']
            if (key, k) not in cache:
                cache[key, k] = load_species(key, k)
            pools[key] = cache[key, k]
        rows = draw_schedule(arm, spec, pools)
        assert not any(r['seqid'] in forbidden or r['seqid'].endswith('M') and '__val' in r['seqid'] for r in rows)
        with open(out / f'schedule_{arm}.jsonl', 'w') as f:
            for r in rows:
                f.write(json.dumps(r) + '\n')
        uniq = {r['id'] for r in rows}
        pos_by_window = {}
        for r in rows:
            if r['id'] not in pos_by_window:
                pos_by_window[r['id']] = positives(pools[r['species']], r['seqid'], r['strand'], r['start'], r['end'])
        exposures = Counter()
        for r in rows:
            for g in pos_by_window[r['id']]:
                exposures[r['species'], g] += 1
        summary[arm] = dict(arm=arm, species=spec['species'], new_species_label_fraction=spec['new_k'], treatment=spec['treatment'], steps=STEPS, batch=BATCH,
            draws=len(rows), unique_windows=len(uniq), draws_per_species=dict(Counter(r['species'] for r in rows)),
            unique_windows_per_species=dict(Counter(i.split('|')[0] for i in uniq)),
            gene_first_fraction=sum(r['mode'] == 'G' for r in rows) / len(rows),
            labeled_genes={k: len(pools[k]['genes']) for k in spec['species']}, eligible_genes={k: pools[k]['eligible'] for k in spec['species']},
            labeled_genes_exposed_unique=len(exposures), labeled_gene_exposures_total=sum(exposures.values()),
            labeled_genes_exposed_unique_per_species=dict(Counter(k[0] for k in exposures)),
            windows_with_positive_fraction=sum(bool(pos_by_window[r['id']]) for r in rows) / len(rows),
            sampling_seed=SAMPLING_SEED, uniform_fraction=UNIFORM_FRACTION)
        print(arm, summary[arm]['unique_windows'], summary[arm]['draws_per_species'], flush=True)
    json.dump(summary, open(sp_path, 'w'), indent=1)


if __name__ == '__main__':
    main()
