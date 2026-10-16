#!/usr/bin/env python3
"""A TRAIN-only metadata census; no DNA, feature, checkpoint or model execution."""
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from src.m28.labels import oriented_parts

ROOT = Path('/home/users/j/jwang/coding-rna')
LABELS = ROOT/'outputs/M28-LABEL-POLICY-R4'
DRAWS = ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl'
ALLOW = {'arabidopsis_thaliana': {'NC_003070.9','NC_003071.7','NC_003076.8'},
         'oryza_sativa': {'NC_089043.1','NC_089044.1','NC_089046.1'}}


def train_records(path):
    # Same pre-parse split filter as the fixed B2 trainer, without importing it.
    with path.open() as f:
        for line in f:
            if '"split": "train"' not in line:
                continue
            row = json.loads(line)
            if row['split'] != 'train':
                raise ValueError('TRAIN metadata filter mismatch')
            yield row


def profile(windows, catalog):
    counts, histogram, exposed = Counter(), Counter(), set()
    maximum, witnesses, max_exons, min_intron = -1, [], 0, None
    for w in windows:
        n = w['valid_bases']
        ids = w['positive_ids']
        local = w['overlapping_transcript_ids']
        if len(ids) != len(set(ids)) or not set(ids) <= set(local):
            raise ValueError('Duplicate or unjoined positive transcript IDs')
        if w['width'] != 24576 or w['end']-w['start'] != n:
            raise ValueError('Frozen window geometry mismatch')
        all_parts = {}
        for ident in local:
            tx = catalog[ident]
            if (tx['seqid'],tx['strand']) != (w['seqid'],w['strand']):
                raise ValueError('TRAIN transcript/window identity mismatch')
            all_parts[ident] = tuple(oriented_parts(tx,w))
        positive = [all_parts[i] for i in ids]
        if len(positive) != len(set(positive)):
            counts['windows_with_duplicate_positive_chain_geometry'] += 1
        for ident, chain in zip(ids,positive):
            tx = catalog[ident]
            if not tx['eligible_positive'] or not tx['training_primary']:
                raise ValueError('Positive transcript eligibility mismatch')
            if not chain or any(not 0 <= a < b <= n for a,b in chain):
                raise ValueError('Positive chain outside oriented valid window')
            if chain[0][0] < 256 or chain[-1][1] > n-256:
                raise ValueError('R4 positive flank mismatch')
            max_exons = max(max_exons,len(chain))
            counts['positive_exon_exposures'] += len(chain)
            counts['positive_chains_over_128_exons'] += int(len(chain)>128)
            for (_,b),(a,_) in zip(chain,chain[1:]):
                gap = a-b
                min_intron = gap if min_intron is None else min(min_intron,gap)
                counts['positive_intron_exposures_below_5bp'] += int(gap<5)
        protected = {p for p in all_parts.values() if p[0][0]>=0 and p[-1][1]<=n}
        other = protected-set(positive)
        unknown = bool(w['chain_unknown_intervals'])
        counts['windows'] += 1
        counts['positive_chain_exposures'] += len(ids)
        counts['windows_with_positive'] += int(bool(ids))
        counts['windows_with_unknown_chain_intervals'] += int(unknown)
        counts['windows_with_full_protected_nonpositive_chains'] += int(bool(other))
        counts['windows_with_either_metadata_protection'] += int(unknown or bool(other))
        counts['reference_empty_metadata_windows'] += int(w['no_annotated_gene_either_strand'])
        counts['zero_positive_windows_with_protection'] += int(not ids and (unknown or bool(other)))
        counts['zero_positive_windows_without_metadata_protection'] += int(not ids and not unknown and not other)
        histogram[len(ids)] += 1
        exposed.update(ids)
        if len(ids)>maximum:
            maximum,witnesses = len(ids),[w['id']]
        elif len(ids)==maximum and len(witnesses)<3 and w['id'] not in witnesses:
            witnesses.append(w['id'])
    return {'counts':dict(counts),'positive_genes_per_window_histogram':{str(k):histogram[k] for k in sorted(histogram)},
            'maximum_positive_genes':maximum,'maximum_witness_window_ids':witnesses,
            'unique_positive_transcript_ids_exposed':len(exposed),
            'maximum_exons_one_positive_chain':max_exons,'minimum_positive_intron_bp':min_intron,
            'metadata_protection_window_fraction':counts['windows_with_either_metadata_protection']/counts['windows']}


def main(out):
    windows,catalogs,by_id = {},{},{}
    for species,allowed in ALLOW.items():
        ws = list(train_records(LABELS/species/'windows_train.jsonl'))
        if not ws or {w['seqid'] for w in ws} != allowed:
            raise ValueError('Frozen TRAIN allowlist mismatch')
        for w in ws:
            if w['species'] != species or w['label_policy'] != 'CDS_assessable_primary_R4' or w['id'] in by_id:
                raise ValueError('TRAIN window grain/policy mismatch')
            by_id[w['id']] = w
        needed = {ident for w in ws for ident in w['overlapping_transcript_ids']}
        catalog = {}
        for tx in train_records(LABELS/species/'transcripts.jsonl'):
            if tx['id'] not in needed:
                continue
            if tx['id'] in catalog or tx['seqid'] not in allowed:
                raise ValueError('TRAIN catalog grain/allowlist mismatch')
            catalog[tx['id']] = tx
        if set(catalog) != needed:
            raise ValueError('TRAIN transcript join coverage mismatch')
        windows[species],catalogs[species] = ws,catalog
    with DRAWS.open() as f:
        draws = [json.loads(line) for line in f if line.strip()]
    if (len(draws)!=1536 or len({d['window_id'] for d in draws})!=1491 or
        Counter(d['species'] for d in draws)!=Counter({s:768 for s in ALLOW})):
        raise ValueError('Original paired draw contract mismatch')
    for d in draws:
        w = by_id[d['window_id']]
        if (w['species'],w['seqid'],w['strand']) != (d['species'],*d['stratum']):
            raise ValueError('Original draw/window identity mismatch')
        if not math.isfinite(d['importance_weight_uniform_within_stratum']) or d['importance_weight_uniform_within_stratum']<=0:
            raise ValueError('Original draw importance mismatch')
    result = {'scope':'TRAIN metadata only; no model performance',
              'label_source':str(LABELS.relative_to(ROOT)),'draw_source':str(DRAWS.relative_to(ROOT)),
              'TRAIN_allowlists':{s:sorted(v) for s,v in ALLOW.items()},
              'views':{},'draw_order_weights_unchanged':True,
              'optimizer_steps':0,'DEV_reference_parsed':False,'DNA_features_checkpoints_read':False,
              'limitations':['metadata protection is a lower bound; callable DNA not inspected',
                             'reference-empty metadata does not certify biological negatives',
                             'query capacity is not selected and DEV/generalization capacity is untested']}
    for species in ALLOW:
        result['views'][species] = {
            'all_natural_TRAIN_windows':profile(windows[species],catalogs[species]),
            'original_768_draws_with_replacement':profile(
                [by_id[d['window_id']] for d in draws if d['species']==species],catalogs[species]),
            'unique_drawn_windows':len({d['window_id'] for d in draws if d['species']==species})}
    q = max(v['all_natural_TRAIN_windows']['maximum_positive_genes'] for v in result['views'].values())
    result['TRAIN_positive_capacity_lower_bound'] = q
    result['teacher_pairs_one_original_pass_if_Q_equals_lower_bound'] = q*sum(
        v['original_768_draws_with_replacement']['counts']['positive_chain_exposures'] for v in result['views'].values())
    out.mkdir(parents=True,exist_ok=True)
    with (out/'summary.json').open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'completed':True,'TRAIN_capacity_lower_bound':q,
                      'natural_windows':{s:len(ws) for s,ws in windows.items()}},allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,required=True)
    main(parser.parse_args().output_dir)
