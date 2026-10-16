#!/usr/bin/env python3
"""A7: fixed paired all-window pooling of completed A6 chain-quality outputs."""
from collections import Counter, defaultdict
import csv
import json
import time
import evaluate_chain_quality as F
from evaluate_pilot import aggregate
from evaluate_window_pool import pool_audit
from src.m28.window_pool import assemble_window_pool
from src.m28.export import as_transcript, write_gff3

Q, A = F.Q, F.A
OUT = Q.ROOT/'outputs/M28-A7-QUALITY-WINDOW-POOL'
PROTOCOL = 'reports/M28-METHOD-RESTART/a7_quality_window_protocol.md'


def require_parent(fitted, decoded, evaluated):
    F.validate_complete(fitted,decoded)
    if (not evaluated['completed'] or evaluated['reference_primary_chains'] != 7728 or
        evaluated['test_setaria_used'] or set(evaluated['methods']) != {'A6'} or
        evaluated['fit'] != fitted or evaluated['protocol'] != Q.PROTOCOL or
        evaluated['inference_receipt']['checkpoint'] != decoded['checkpoint']):
        raise ValueError('Completed evaluation of the exact fixed A6 parent required')


def assemble_scope(path,scope):
    assembled = assemble_window_pool(F.records(path,scope),scope['windows'],scope['strand'])
    if any(assembled['counts'][k] != scope['assembly_counts'][k]
           for k in ('free_records','free_unique')):
        raise ValueError('A7 must retain exactly the completed A6 free pool')
    return assembled


def paired_references(refs,before,after):
    rows = []
    for species,references in refs.items():
        for ref in references:
            key = A.E.S.key(ref)
            rows.append({'species':species,'id':ref['id'],'gene_id':ref['gene_id'],
                         'seqid':ref['seqid'],'strand':ref['strand'],
                         'CDS_intervals':json.dumps([[a,b] for a,b,_ in ref['CDS']],separators=(',',':')),
                         'unsupported_strata':';'.join(ref.get('M28_unsupported_strata_not_excluded',[])),
                         'A6_selected':key in before,'A7_selected':key in after})
    return rows


def membership(rows):
    counts = Counter()
    for r in rows:
        a,b = r['A6_selected'],r['A7_selected']
        counts['both' if a and b else 'A7_only' if b else 'A6_only' if a else 'neither'] += 1
    return {k:counts[k] for k in ('both','A7_only','A6_only','neither')}


def main():
    began = time.perf_counter()
    parent_eval = Q.ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL'
    if OUT.exists():
        raise FileExistsError('Do not overwrite an earlier A7 run')
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (Q.FIT,Q.INFER,parent_eval)):
        raise ValueError('Terminal A6 fit/inference/evaluation required')
    fitted,decoded,parent = [json.loads((p/'summary.json').read_text()) for p in (Q.FIT,Q.INFER,parent_eval)]
    require_parent(fitted,decoded,parent)
    OUT.mkdir()
    lengths = {q:n for qs in Q.I.DEV_LENGTHS.values() for q,n in qs.items()}
    predictions,scopes = defaultdict(list),{}
    # Pooling and GFF export are completed before any reference is opened.
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        assembled = assemble_scope(Q.INFER/scope['file'],scope)
        scopes[species+'|'+strand] = (scope,assembled)
        for i,row in enumerate(assembled['predictions']):
            predictions[species].append(as_transcript(row['chain'],seqid,strand,
                'M28_A7_'+str(scope['ordinal'])+'_'+strand+'_'+str(i),row['score']))
    txs = [t for p in predictions.values() for t in p]
    gff = OUT/'A7.gff3'
    write_gff3(gff,txs,lengths,'M28_A7')
    saved = A.E.S.E.primary_transcripts(A.E.S.E.parse_annotation(gff,lengths),complete_only=False)
    A.validate_gff_roundtrip(saved,txs)
    before_txs = A.E.S.E.primary_transcripts(
        A.E.S.E.parse_annotation(Q.INFER/'A6.gff3',lengths),complete_only=False)
    refs,isoforms,parents,annotations,sequences,_ = A.reference_context()
    stages = {}
    for key,(scope,assembled) in scopes.items():
        species,strand = scope['species'],scope['strand']
        audit = pool_audit(assembled,[r for r in refs[species] if r['strand']==strand])
        original = parent['methods']['A6']['stages'][key]
        if audit['exact_reference_hits']['free_unique'] != original['exact_reference_hits']['free_unique']:
            raise ValueError('A7 reference free coverage differs from A6 source')
        audit.update(species=species,seqid=scope['seqid'],strand=strand)
        audit['A6_owner_control'] = original
        stages[key] = audit
    before,after = {A.E.S.key(t) for t in before_txs},{A.E.S.key(t) for t in txs}
    paired = paired_references(refs,before,after)
    if len(paired) != 7728:
        raise ValueError('Do not drop primary references')
    current = {}
    for species,seqid,length,_ in A.E.S.SCOPE:
        part = [r for r in paired if r['species']==species]
        baseline = parent['methods']['A6']['per_species'][species]['exact_chain']
        baseline_txs = [t for t in before_txs if t['seqid']==seqid]
        if (sum(r['A6_selected'] for r in part) != baseline['tp'] or
            len({A.E.S.key(t) for t in baseline_txs}) != baseline['predicted']):
            raise ValueError('A6 parent GFF and completed metric disagree')
        measured = F.score_species(predictions[species],refs[species],isoforms[species],
                                   length,annotations[species],parents[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        failures = Counter()
        for t in predictions[species]:
            failures.update(Q.I.M.chain_checks(t,sequences[species][seqid])[0])
        measured['CDS_structure_failure_counts_nonexclusive'] = dict(failures)
        if measured['exact_chain']['tp'] != sum(r['A7_selected'] for r in part):
            raise ValueError('Paired reference membership and A7 TP disagree')
        current[species] = measured
    pooled = aggregate(current)
    if sum(s['exact_reference_hits']['selected_unique'] for s in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('A7 stage TP and final TP disagree')
    report = {'completed':True,'experiment':'M28-A7-QUALITY-WINDOW-POOL','protocol':PROTOCOL,
              'status':'complete_paired_DEV_decoder_readout_not_paper_or_generalization_PASS',
              'reference_primary_chains':7728,'reference_filters_added':False,'test_setaria_used':False,
              'reference_used_in_candidate_pooling':False,'new_training':False,'new_model_inference':False,
              'optimizer_steps_this_run':0,'tau':0.,
              'pool_policy':'max_unchanged_A6_chain_quality_gain_across_all_generating_windows',
              'selector':'unchanged_positive_gain_same_strand_span_nonoverlap',
              'parent_checkpoint':fitted['final_checkpoint'],'source_inference':str(Q.INFER.relative_to(Q.ROOT)),
              'source_evaluation':str(parent_eval.relative_to(Q.ROOT)),
              'methods':{'A7':{'per_species':current,'pooled':pooled,'stages':stages,
                              'GFF3':str(gff.relative_to(Q.ROOT))}},
              'A6_parent_control':parent['methods']['A6'],
              'paired_reference_membership':{'pooled':membership(paired),
                    'per_species':{s:membership([r for r in paired if r['species']==s]) for s in refs}},
              'predeclared_matched_owner_gate':A.continuation_gate(current,parent['methods']['A6']['per_species']),
              'A6_original_gate_unchanged':parent['predeclared_pilot_gate'],
              'limitations':'unchanged A6 generated pool/quality score, only ownership-policy contrast; cannot create absent chains; max over windows is not a calibrated biological probability',
              'seconds':time.perf_counter()-began}
    with (OUT/'references.tsv').open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(paired[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(paired)
    Q.R.write_summary(OUT,report)


if __name__=='__main__':
    main()
