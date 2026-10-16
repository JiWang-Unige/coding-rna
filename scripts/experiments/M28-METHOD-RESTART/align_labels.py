#!/usr/bin/env python3
"""Align full-chain targets to the frozen ruler, retaining original sampling law."""
import argparse
import json
from collections import Counter
from pathlib import Path
import common_ruler as E
import manifest as M

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def relabel_window(window,catalog):
    a,b=window['start'],window['end']
    local=[catalog[i] for i in window['overlapping_transcript_ids']]
    fully=[t for t in local if a<=t['span'][0] and t['span'][1]<=b]
    pos=[t['id'] for t in fully if t['eligible_positive'] and a+M.FLANK<=t['span'][0] and t['span'][1]<=b-M.FLANK]
    unknown,fine=[],[]
    for t in local:
        interval=M.orient_interval(max(a,t['span'][0]),min(b,t['span'][1]),a,b,window['strand'])
        if not t['CDS_complete_supported'] or t['overlapping_primary']:
            fine.append(interval);unknown.append(interval)
        elif t not in fully:
            unknown.append(interval)
    return dict(window,positive_ids=pos,chain_unknown_intervals=M.merged(unknown),
                fine_label_unknown_intervals=M.merged(fine),label_policy='CDS_assessable_primary_R4')

def main(species,out):
    out.mkdir(parents=True,exist_ok=True)
    meta=M.C.metadata(M.C.yaml.safe_load(M.C.CONFIG.read_text()))[species]
    seqs=M.read_allowed(meta['path']/'genome.fa',meta['lengths'])
    coding=E.S.E.parse_annotation(meta['path']/'reference.gff3',meta['lengths'],protein_coding_only=True)
    raw_cds={}
    for seqid in meta['lengths']:
        with (meta['path']/'reference.gff3').open() as f:
            raw_cds.update(E.P.raw_features(f,seqid)[2])
    refs,_=E.select_reference(coding,raw_cds,seqs)
    if any(not E.R7.phase_chain_consistent(t) for t in refs):
        raise ValueError('Assessable reference phase mismatch')
    reference_by_gene={t['gene_id']:t for t in refs}
    local_primary={t['gene_id']:t for t in E.S.E.primary_transcripts(coding,complete_only=False)}
    changed=[g for g,t in reference_by_gene.items() if E.S.key(t)!=E.S.key(local_primary[g])]
    # Genes with no assessable complete isoform retain certain local evidence only.
    local_primary.update(reference_by_gene)
    primary_ids={t['id'] for t in local_primary.values()}
    reference_ids={t['id'] for t in refs}
    overlaps=set()
    for seqid in meta['seqids']:
        overlaps |= M.C.overlapping_ids([t for t in local_primary.values() if t['seqid']==seqid])
    source=M.C.ROOT/'outputs/M28-DATA-MANIFEST-R2'/species
    catalog={t['id']:t for t in rows(source/'transcripts.jsonl')}
    for t in catalog.values():
        t['training_primary']=t['id'] in primary_ids
        t['CDS_assessable_primary']=t['id'] in reference_ids
        t['local_only_fallback_primary']=t['training_primary'] and t['gene_id'] not in reference_by_gene
        t['overlapping_primary']=t['id'] in overlaps
        t['eligible_positive']=t['CDS_assessable_primary'] and t['CDS_complete_supported'] and not t['overlapping_primary']
    with (out/'transcripts.jsonl').open('x') as f:
        for t in catalog.values(): f.write(json.dumps(t)+'\n')
    windows={}
    split_summary={}
    for split in ('train','val'):
        new=[relabel_window(w,catalog) for w in rows(source/('windows_'+split+'.jsonl'))]
        with (out/('windows_'+split+'.jsonl')).open('x') as f:
            for w in new: f.write(json.dumps(w)+'\n')
        windows.update({w['id']:w for w in new})
        split_summary[split]={'windows':len(new),'CDS_assessable_primaries':sum(t['split']==split and t['CDS_assessable_primary'] for t in catalog.values()),
                              'eligible_positive':sum(t['split']==split and t['eligible_positive'] for t in catalog.values()),
                              'local_only_fallback_primaries':sum(t['split']==split and t['local_only_fallback_primary'] for t in catalog.values()),
                              'primary_chain_changes':sum(catalog[reference_by_gene[g]['id']]['split']==split for g in changed)}
    frozen=rows(M.C.ROOT/'outputs/M28-COMMON-RULER-R3/result'/(species+'.references.jsonl'))
    assert {t['id'] for t in frozen}=={t['id'] for t in catalog.values() if t['split']=='val' and t['CDS_assessable_primary']}
    draws=[d for d in rows(M.C.ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl') if d['species']==species]
    assert len(draws)==768
    # Do not resample or recompute q from the new labels: original q is the actual draw law.
    exposure=Counter(i for d in draws for i in windows[d['window_id']]['positive_ids'])
    old_windows={w['id']:w for w in rows(source/'windows_train.jsonl')}
    changed_draws=sum(windows[d['window_id']]['positive_ids']!=old_windows[d['window_id']]['positive_ids'] for d in draws)
    report={'species':species,'label_policy':'longest_CDS_assessable_primary_then_training_support_mask',
            'local_only_fallback':'longest_all when no assessable isoform; certain local/junction evidence only, never complete-chain positive',
            'draw_source':'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl',
            'draw_order_q_weights_unchanged':True,'draws':len(draws),
            'positive_list_changed_draws':changed_draws,'unique_positive_chains_exposed':len(exposure),
            'primary_change_gene_ids':changed,'splits':split_summary,
            'DEV_reference_ID_set_matches_frozen_ruler':True,'optimizer_steps':0,'test_setaria_used':False}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--species',choices=('arabidopsis_thaliana','oryza_sativa'),required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.species,a.output_dir)
