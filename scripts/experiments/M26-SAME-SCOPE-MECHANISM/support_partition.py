#!/usr/bin/env python3
"""Reference-assisted membership diagnostic of the unchanged R1 candidate graph."""
import argparse
import csv
import json
import sys
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import same_scope as S

sys.path.insert(0, str(S.ROOT/'scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION'))
import run_graph_core as G


def witness(options):
    previous, chosen = -1, []
    for values in options:
        values = sorted(set(values))
        index = bisect_right(values, previous)
        if index == len(values):
            return None
        previous = values[index]
        chosen.append(previous)
    return chosen


def membership(chain, index, orf):
    if sum(b-a for a,b in chain) < 6 or any(b+4>c for (_,b),(c,_) in zip(chain,chain[1:])):
        return False
    options, fragments = [], []
    for i,(a,b) in enumerate(chain):
        left = 'start' if i == 0 else 'acceptor'
        right = 'stop' if i == len(chain)-1 else 'donor'
        options.append(index.get((a,b,left,right), []))
        fragments.append((a,b,left,right))
    chosen = witness(options)
    if chosen is None:
        return False
    state = ''
    for rid,(a,b,left,right) in zip(chosen,fragments):
        state = orf.advance(state, G.Exon(rid,a,b,left,right,0.0))
        if state is None:
            return False
    return state == ''


def simultaneous_bound(txs):
    groups = defaultdict(list)
    for t in txs:
        groups[(t['seqid'],t['strand'])].append((t['CDS'][0][0],t['CDS'][-1][1]))
    total = 0
    for spans in groups.values():
        last = -1
        for a,b in sorted(spans,key=lambda span:(span[1],span[0])):
            if a >= last:
                total += 1
                last = b
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'support_partition.json').exists():
        raise FileExistsError('Do not overwrite a completed diagnostic')
    with (S.ROOT/'outputs/M25R-E1-BLOCK-ASSEMBLY-REACHABILITY/references.jsonl').open() as f:
        oracle = {(r['species'],r['seqid'],r['strand'],r['transcript_id']):r for r in map(json.loads,f)}
    assert len(oracle) == 6450
    ledger, references, prediction_keys = [], [], defaultdict(set)
    for species,seqid,length,expected in S.SCOPE:
        seq = S.V.read_fasta(S.ROOT/'data/m1_screen'/species/'genome.fa')[seqid]
        assert len(seq) == length
        ann = S.E.parse_annotation(S.ROOT/'data/m1_screen'/species/'reference.gff3', {seqid:length}, protein_coding_only=True)
        refs = S.E.primary_transcripts(ann)
        assert len(refs) == expected
        references.extend(refs)
        for method,template in S.METHODS.items():
            path = S.ROOT/template.format(species=species)
            txs = S.E.primary_transcripts(S.E.parse_annotation(path,{seqid:length}),complete_only=False)
            prediction_keys[method].update(S.key(t) for t in txs)
        for strand in ('+','-'):
            orientation = seq if strand == '+' else S.M.reverse_complement(seq)
            stem = f'{seqid}_{"plus" if strand == "+" else "minus"}'
            cache = S.ROOT/'outputs/M25R-E1-RUN-GRAPH-R1/scores'
            region = np.load(cache/f'{stem}_region.npy',mmap_mode='r')
            boundary = np.load(cache/f'{stem}_boundary.npy',mmap_mode='r')
            assert region.shape == (length,3) and boundary.shape == (length,4)
            assert region.dtype == boundary.dtype == np.float16
            exons,prefix,counts = G.candidates(orientation,region,boundary)
            index = defaultdict(list)
            for e in exons:
                index[(e.a,e.b,e.left,e.right)].append(e.run)
            orf = G.ORF(orientation)
            for t in refs:
                if t['strand'] != strand:
                    continue
                coords = [(a,b) for a,b,_p in t['CDS']]
                oriented = coords if strand == '+' else [(length-b,length-a) for a,b in reversed(coords)]
                old = oracle[(species,seqid,strand,t['id'])]
                o0 = old['structurally_valid_chain_in_graph']
                o1 = membership(oriented,index,orf)
                exact = S.key(t) in prediction_keys['M25R_B']
                assert not o1 or o0, 'Actual support outside relaxed oracle; stop attribution'
                assert not exact or o1, 'B exact chain declared unreachable; stop attribution'
                category = ('B_exact' if exact else 'actual_graph_not_selected' if o1 else
                            'actual_boundary_filter_loss' if o0 else 'outside_relaxed_graph')
                ledger.append({'species':species,'seqid':seqid,'strand':strand,'transcript':t['id'],
                               'category':category,'O0':bool(o0),'O1':bool(o1),
                               'canonical':old['canonical_ORF'],'original_category':old['original_category'],
                               'A_exact':S.key(t) in prediction_keys['M25R_A'],
                               **{method+'_exact':S.key(t) in prediction_keys[method]
                                  for method in ('M25R_B','ANNEVO','Helixer','Tiberius')}})
            print(species,strand,counts,flush=True)
            del exons,prefix,index,orf,region,boundary
    assert len(ledger) == 6450 and sum(r['M25R_B_exact'] for r in ledger) == 3351
    assert sum(r['A_exact'] for r in ledger) == 1389 and sum(r['O0'] for r in ledger) == 3585
    actual_ids = {(r['seqid'],r['strand'],r['transcript']) for r in ledger if r['O1']}
    actual_refs = [t for t in references if (t['seqid'],t['strand'],t['id']) in actual_ids]
    summary = {'status':'fixed_graph_reference_assisted_development_diagnostic','setaria_access':False,
               'reference_total':len(ledger),'partition':dict(Counter(r['category'] for r in ledger)),
               'O1_individually_reachable':len(actual_refs),
               'O1_max_nonoverlapping_reference_subset':simultaneous_bound(actual_refs),
               'A_B_truth_transitions':dict(Counter(f"A{int(r['A_exact'])}_B{int(r['M25R_B_exact'])}" for r in ledger)),
               'per_species':{},'baseline_recovery_by_partition':{},'by_original_category':{}}
    for species,_q,_n,_c in S.SCOPE:
        summary['per_species'][species] = dict(Counter(r['category'] for r in ledger if r['species']==species))
    for category in sorted({r['category'] for r in ledger}):
        rows = [r for r in ledger if r['category']==category]
        summary['baseline_recovery_by_partition'][category] = {'reference_count':len(rows),
            **{m:sum(r[m+'_exact'] for r in rows) for m in ('ANNEVO','Helixer','Tiberius')}}
    for original in sorted({r['original_category'] for r in ledger}):
        summary['by_original_category'][original] = dict(Counter(r['category'] for r in ledger if r['original_category']==original))
    with (out/'support_partition.json').open('x') as f:
        json.dump(summary,f,indent=2)
        f.write('\n')
    with (out/'support_ledger.tsv').open('x') as f:
        writer = csv.DictWriter(f,fieldnames=list(ledger[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(ledger)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
