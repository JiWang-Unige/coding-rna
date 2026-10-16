#!/usr/bin/env python3
"""One R1 execution, split into reference-free prediction and evaluation processes.

infer -> replay -> evaluate A -> decode -> evaluate B.
No fitting, held-out access, or reference input to graph decoding.
CPU recovery may use an independent output directory and existing saved scores.
"""
import argparse
from collections import Counter,defaultdict
import json
import time
from pathlib import Path
import numpy as np
import yaml
import redecode_error_decomposition as D
import run_graph_core as G

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'outputs/M25R-E1-RUN-GRAPH-R1'
SCORE_SOURCE=None
CONFIG=ROOT/'configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml'
CHECKPOINT=ROOT/'outputs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/checkpoints/epoch_1.pt'
FROZEN=ROOT/'outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4/epoch_1'
SCOPE=(('arabidopsis_thaliana','NC_003074.8',23459830),('oryza_sativa','NC_089041.1',29936421))
THRESHOLDS={'region':0.4,'start':0.5,'stop':0.5,'donor':0.1,'acceptor':0.1}


def read(path):
    return json.loads(path.read_text())


def save(name,payload):
    with (OUT/name).open('x') as f:
        json.dump(payload,f,indent=2)
        f.write('\n')


def sequences():
    result={}
    for species,seqid,length in SCOPE:
        seq=D.screen_data.read_fasta(str(ROOT/'data/m1_screen'/species/'genome.fa'))[seqid].upper()
        if len(seq)!=length:
            raise AssertionError('approved validation chromosome length changed')
        result[seqid]=seq
    return result


def score_path(seqid,strand,name,run_dir=None):
    return (run_dir or OUT)/'scores'/f'{seqid}_{"plus" if strand=="+" else "minus"}_{name}.npy'


def load_scores(seqid,strand):
    arrays=[np.load(score_path(seqid,strand,name,SCORE_SOURCE),mmap_mode='r')
            for name in ('region','boundary','phase')]
    length=next(n for _s,q,n in SCOPE if q==seqid)
    if any(a.shape!=(length,w) or a.dtype!=np.float16 for a,w in zip(arrays,(3,4,4))):
        raise AssertionError('full saved logits have wrong shape or precision')
    return arrays


def infer():
    import torch,transformers,peft
    if (OUT/'STATUS').read_text().strip()!='RUNNING_INFERENCE':
        raise AssertionError('inference stage entered outside approved state')
    if not (OUT/'synthetic_tests.xml').is_file():
        raise AssertionError('synthetic validation evidence missing')
    config=yaml.safe_load(CONFIG.read_text())
    if (torch.__version__,transformers.__version__,peft.__version__)!=('2.5.1','4.49.0','0.19.1'):
        raise AssertionError('frozen inference environment changed')
    if config['seed']!=0 or config['model']['window_bp']!=6144:
        raise AssertionError('frozen config changed')
    seqs=sequences()
    expected=sum(2*((n+6143)//6144) for _s,_q,n in SCOPE)
    assert expected==17384
    save('resolved_run.json',{'contract':'M25R-E1-RUN-GRAPH-R1','checkpoint':str(CHECKPOINT),
         'scope':SCOPE,'windows':expected,'thresholds':THRESHOLDS,'epoch':1,'enumeration_order':601,
         'training':False,'setaria_access':False,'reference_in_inference_or_decoder':False,
         'saved_dtype':'float16','score_accumulation':'float64','weights':'all_one'})
    (OUT/'scores').mkdir()
    tokenizer,backbone,heads,forward,device,k=D.load_inference_model(config)
    D.load_checkpoint(CHECKPOINT,backbone,heads)
    for module in (backbone,heads):
        module.eval()
        for p in module.parameters(): p.requires_grad_(False)
    started=time.monotonic()
    count=0
    with torch.no_grad(),(OUT/'window_runs.jsonl').open('x') as log:
        for species,seqid,length in SCOPE:
            for strand in ('+','-'):
                sequence=seqs[seqid] if strand=='+' else D.m25.reverse_complement(seqs[seqid])
                arrays=[np.lib.format.open_memmap(score_path(seqid,strand,name),mode='w+',dtype=np.float16,
                        shape=(length,width)) for name,width in (('region',3),('boundary',4),('phase',4))]
                for start in range(0,length,6144):
                    end=min(start+6144,length)
                    padded=sequence[start:end]+'A'*(6144-(end-start))
                    ids=D._tokenize_window(tokenizer,D._clean(padded),6144,k).unsqueeze(0).to(device)
                    attention=torch.ones_like(ids)
                    nucleotide=torch.from_numpy(D.m25._one_hot(padded)).unsqueeze(0).to(device)
                    with torch.autocast(device_type='cuda',dtype=torch.bfloat16):
                        values=forward(ids,attention,nucleotide)
                    for target,value in zip(arrays,values):
                        value=value[0,:end-start].float().cpu().numpy().astype(np.float16)
                        if not np.isfinite(value).all(): raise AssertionError('nonfinite inference logits')
                        target[start:end]=value
                    count+=1
                    elapsed=time.monotonic()-started
                    log.write(json.dumps({'window':count,'species':species,'seqid':seqid,'strand':strand,
                                         'start':start,'end':end,'elapsed_seconds':elapsed})+'\n')
                    if count%50==0:
                        log.flush()
                        print(f'window={count}/{expected} elapsed={elapsed:.1f}s',flush=True)
                for a in arrays: a.flush()
                del arrays
                print(f'orientation_saved={seqid}/{strand}',flush=True)
    assert count==expected
    save('inference_complete.json',{'windows_forwarded':count,'elapsed_seconds':time.monotonic()-started,
                                  'cuda_max_allocated_bytes':torch.cuda.max_memory_allocated()})


def export(arm,predictions):
    payload=[{'species':species,'seqid':seqid,'models':models} for (species,seqid),models in predictions.items()]
    save(f'{arm}_models.json',payload)
    if arm=='A':
        D.m25._write_gff3(str(OUT/'A_predictions.gff3'),{q:m for (_s,q),m in predictions.items()},'M25R_R1_A')
        return
    # Codons can span CDS joins. Emit their actual CDS pieces, not a contiguous
    # genomic 3bp feature that incorrectly includes intronic bases.
    serial=0
    with (OUT/'B_predictions.gff3').open('x') as f:
        f.write('##gff-version 3\n')
        for (_species,seqid),models in sorted(predictions.items()):
            for m in sorted(models,key=lambda m:(m['cds'][0][0],m['strand'],m['cds'])):
                serial+=1
                gene=f'M25R_R1_gene_{serial:07d}'
                tx=gene+'.t1'
                def line(feature,a,b,phase,attrs):
                    f.write(f'{seqid}\tM25R_R1_B\t{feature}\t{a+1}\t{b}\t.\t{m["strand"]}\t{phase}\t{attrs}\n')
                line('gene',m['cds'][0][0],m['cds'][-1][1],'.',f'ID={gene}')
                line('mRNA',m['cds'][0][0],m['cds'][-1][1],'.',f'ID={tx};Parent={gene}')
                for (a,b),p in zip(m['cds'],m['phase']): line('CDS',a,b,p,f'Parent={tx}')
                for feature,parts in codon_pieces(m['cds'],m['strand']).items():
                    for a,b in parts: line(feature,a,b,0,f'Parent={tx}')


def codon_pieces(chain,strand):
    ordered=list(chain) if strand=='+' else list(reversed(chain))
    def take(parts,from_start):
        remaining,result=3,[]
        for a,b in parts:
            width=min(remaining,b-a)
            result.append((a,a+width) if from_start else (b-width,b))
            remaining-=width
            if remaining==0: break
        assert remaining==0
        return sorted(result)
    return {'start_codon':take(ordered,strand=='+'),
            'stop_codon':take(list(reversed(ordered)),strand=='-')}


def replay():
    assert read(OUT/'inference_complete.json')['windows_forwarded']==17384
    seqs=sequences()
    predictions={}
    for species,seqid,length in SCOPE:
        models=[]
        for strand in ('+','-'):
            seq=seqs[seqid] if strand=='+' else D.m25.reverse_complement(seqs[seqid])
            arrays=load_scores(seqid,strand)
            _states,_traces,_prefilter,emitted=D.trace_decode_orientation(seq,*arrays,THRESHOLDS,
                                                                         reverse_mapped=strand=='-')
            models.extend(D.m25._map_model_from_orientation(m,length,strand) for m in emitted)
            print(f'A_replay={seqid}/{strand} emitted={len(emitted)}',flush=True)
        predictions[(species,seqid)]=models
    export('A',predictions)
    save('A_prediction_complete.json',{'n_models':sum(map(len,predictions.values()))})


def decode():
    assert read(OUT/'A_replay_check.json')['passed'] is True
    seqs=sequences()
    predictions={}
    statistics=[]
    for species,seqid,length in SCOPE:
        models=[]
        for strand in ('+','-'):
            sequence=seqs[seqid] if strand=='+' else D.m25.reverse_complement(seqs[seqid])
            region,boundary,_phase=load_scores(seqid,strand)
            exons,prefix,stats=G.candidates(sequence,region,boundary)
            path=G.solve(sequence,exons,prefix,strand)
            oriented=G.models(path)
            models.extend(D.m25._map_model_from_orientation(m,length,strand) for m in oriented)
            stats.update(seqid=seqid,strand=strand,objective=path.score,total_gene_span=path.span,
                         genes=len(oriented))
            statistics.append(stats)
            print(json.dumps(stats),flush=True)
            del exons,prefix,path,region,boundary,_phase
        predictions[(species,seqid)]=models
    export('B',predictions)
    save('B_prediction_complete.json',{'orientations':statistics,'reference_files_read':False,
         'n_models':sum(map(len,predictions.values())),'phase_gate':False,
         'ORF_suffix_states':21,'initial_ATG_partial_states':2})


def audited_transcripts(path,seqs,expected_models,split_codons=False):
    lengths={q:len(s) for q,s in seqs.items()}
    if not split_codons:
        audit=D.audit_gff3(path,seqs,lengths,expected_models)
    else:
        # Retain original sequence/phase/splice checks, but derive codon feature
        # coordinates independently from concatenated CDS coordinates.
        genes,txs,children={}, {}, defaultdict(list)
        for number,line in enumerate(path.open(),1):
            if line.startswith('#') or not line.strip(): continue
            q,_src,feature,a,b,_score,strand,phase,attrs=line.rstrip().split('\t')
            a,b=int(a)-1,int(b)
            if q not in lengths or not 0<=a<b<=lengths[q] or strand not in {'+','-'}:
                raise AssertionError(f'invalid output coordinates at line {number}')
            attrs=D.m25_eval.m24.parse_attrs(attrs)
            if feature=='gene': genes[attrs['ID']]=(q,a,b,strand)
            elif feature=='mRNA': txs[attrs['ID']]=(attrs['Parent'],q,a,b,strand)
            elif feature in {'CDS','start_codon','stop_codon'}: children[attrs['Parent']].append((feature,a,b))
        if set(children)!=set(txs) or len(genes)!=len(txs) or len(txs)!=expected_models:
            raise AssertionError('gene/transcript/child linkage count mismatch')
        for tx,(parent,q,a,b,strand) in txs.items():
            if genes.get(parent)!=(q,a,b,strand): raise AssertionError('gene parent mismatch')
        annotation=D.parse_annotation(str(path),lengths,protein_coding_only=False)
        parsed=D.primary_transcripts(annotation)
        if len(parsed)!=expected_models: raise AssertionError('primary parsed output count mismatch')
        codons=D.m25_eval.codon_features(path)
        ledger=[]
        failures=Counter()
        for tx in parsed:
            components={'parent_linkage':True,**D.m25_eval.structural_validity_components(tx,seqs,codons)}
            # Independent per-base coordinates for the first and last 3 CDS bp.
            coords=[]
            for a,b,_p in D.m25_eval.ordered_cds(tx):
                coords.extend(range(a,b) if tx['strand']=='+' else range(b-1,a-1,-1))
            for feature,wanted in (('start_codon',coords[:3]),('stop_codon',coords[-3:])):
                observed=[]
                for q,a,b,strand in codons.get(tx['id'],{}).get(feature,[]):
                    if q!=tx['seqid'] or strand!=tx['strand']: raise AssertionError('codon strand/seqid mismatch')
                    observed.extend(range(a,b))
                components[feature+'_feature']=sorted(observed)==sorted(wanted) and len(observed)==3
            failed=[k for k,v in components.items() if not v]
            failures.update(failed)
            ledger.append({'transcript':tx['id'],'valid':not failed,'components':components,'failed_components':failed})
        valid=sum(r['valid'] for r in ledger)
        audit={'n_emitted':expected_models,'n_checked':len(parsed),'audit_coverage':1.0,
               'valid_transcripts':valid,'invalid_transcripts':expected_models-valid,
               'validity_fraction':valid/expected_models if expected_models else 'not_applicable',
               'component_failure_counts':dict(failures),'transcript_ledger':ledger,
               'codon_feature_rule':'first/last 3 spliced CDS bases, split records permitted'}
    save(path.stem+'_structural_audit.json',audit)
    if audit['invalid_transcripts'] or audit['audit_coverage']!=1.0:
        raise AssertionError('independent structural audit failed; stop R1')
    return audit


def chainset(path,lengths):
    annotation=D.parse_annotation(str(path),lengths,protein_coding_only=False)
    return {(tx['seqid'],tx['strand'],tuple((a,b,p) for a,b,p in tx['CDS']))
            for tx in D.primary_transcripts(annotation)}


def load_predictions(path):
    rows=read(path)
    predictions={(r['species'],r['seqid']):r['models'] for r in rows}
    # JSON erases tuple types. The unchanged frozen metric stores chains in sets.
    for models in predictions.values():
        for model in models:
            model['cds']=[tuple(interval) for interval in model['cds']]
    return predictions


def evaluate(arm):
    # This is the only stage that opens development reference GFFs. Its caller
    # starts it only after prediction has exited and completely saved its GFF.
    assert (OUT/f'{arm}_prediction_complete.json').is_file()
    predictions=load_predictions(OUT/f'{arm}_models.json')
    seqs=sequences()
    lengths={q:len(s) for q,s in seqs.items()}
    path=OUT/f'{arm}_predictions.gff3'
    audit=audited_transcripts(path,seqs,sum(map(len,predictions.values())),split_codons=arm=='B')
    exported=chainset(path,lengths)
    in_memory={(q,m['strand'],tuple((a,b,str(p)) for (a,b),p in zip(m['cds'],m['phase'])))
               for (_s,q),ms in predictions.items() for m in ms}
    if exported!=in_memory or len(exported)!=sum(map(len,predictions.values())):
        raise AssertionError('saved GFF differs from decoded chains or has duplicates')
    config=yaml.safe_load(CONFIG.read_text())
    species=D.load_species(ROOT,config)
    refs,ref_lengths=D.validation_truth(species)
    if set(refs)!={(s,q) for s,q,_n in SCOPE} or sum(map(len,refs.values()))!=6450:
        raise AssertionError('frozen complete-primary reference scope changed')
    metrics=D.m25._validation_metrics(predictions,refs,ref_lengths)
    metrics['structurally_valid_complete_fraction']=audit['validity_fraction'] if sum(map(len,predictions.values())) else 0.0
    if arm=='A':
        old=read(FROZEN/'diagnostic.json')
        if old['frozen_tuple']!={**THRESHOLDS,'epoch':1,'enumeration_order':601}:
            raise AssertionError('R4 frozen tuple changed')
        same=exported==chainset(FROZEN/'replayed_predictions.gff3',lengths)
        repro={name:{'replayed':metrics[name],'frozen':r['replayed'],
                     'absolute_error':abs(metrics[name]-r['replayed'])} for name,r in old['reproduction'].items()}
        passed=same and all(r['absolute_error']<=1e-5 for r in repro.values())
        save('A_replay_check.json',{'passed':passed,'full_chain_and_phase_set_identical':same,
                                  'metrics':repro,'n_models':len(exported),'audit_invalid':audit['invalid_transcripts']})
        if not passed: raise AssertionError('A replay mismatch; B must not run')
        print('A_FULL_REPLAY_PASS',flush=True)
        (OUT/'STATUS').write_text('A_FULL_REPLAY_PASS\n')
        return
    # No tuning follows these real metrics. Necessary single-model gates only.
    errors=D.prediction_error_summary(refs,D.model_transcripts(predictions))
    per_species={}
    for name in sorted(species):
        ps={k:v for k,v in predictions.items() if k[0]==name}
        rs={k:v for k,v in refs.items() if k[0]==name}
        ls={k:v for k,v in ref_lengths.items() if k[0]==name}
        per_species[name]={'metrics':D.m25._validation_metrics(ps,rs,ls),
                           'prediction_errors':D.prediction_error_summary(rs,D.model_transcripts(ps))}
    gate=config['success_gate']
    checks={
        'exact_CDS_interval_F1':metrics['exact_CDS_interval_F1']>=gate['exact_CDS_interval_F1_min'],
        'exact_CDS_chain_F1':metrics['exact_CDS_chain_F1']>=gate['exact_CDS_chain_F1_min'],
        'exact_coding_gene_F1':metrics['exact_CDS_chain_F1']>=gate['exact_coding_gene_F1_min'],
        'matched_gene_strand_accuracy':isinstance(errors['matched_gene_strand']['accuracy'],(float,int)) and errors['matched_gene_strand']['accuracy']>=gate['matched_gene_strand_accuracy_min'],
        'exact_matched_CDS_phase_accuracy':isinstance(errors['exact_matched_CDS_phase']['accuracy'],(float,int)) and errors['exact_matched_CDS_phase']['accuracy']>=gate['exact_matched_CDS_phase_accuracy_min'],
        'structural_validity':metrics['structurally_valid_complete_fraction']>=gate['structurally_valid_complete_transcript_fraction_min'],
        'intergenic_FPR':metrics['intergenic_FPR']<=gate['intergenic_FPR_max'],
        'gene_count_ratio':gate['gene_count_ratio'][0]<=metrics['gene_count_ratio']<=gate['gene_count_ratio'][1]}
    stop=config['immediate_stop_gate']
    immediate={
        'exact_CDS_interval_F1':metrics['exact_CDS_interval_F1']<stop['exact_CDS_interval_F1_below'],
        'exact_CDS_chain_F1':metrics['exact_CDS_chain_F1']<stop['exact_CDS_chain_F1_below'],
        'exact_coding_gene_F1':metrics['exact_CDS_chain_F1']<stop['exact_coding_gene_F1_below'],
        'intergenic_FPR':metrics['intergenic_FPR']>stop['intergenic_FPR_above'],
        'gene_count_ratio':not stop['gene_count_ratio_outside'][0]<=metrics['gene_count_ratio']<=stop['gene_count_ratio_outside'][1],
        'structural_validity':metrics['structurally_valid_complete_fraction']<stop['structurally_valid_complete_transcript_fraction_below']}
    a=chainset(OUT/'A_predictions.gff3',lengths)
    save('A_to_B_chain_changes.json',{'unchanged':len(a&exported),'added':sorted(exported-a),'removed':sorted(a-exported)})
    save('B_evaluation.json',{'metrics':metrics,'prediction_errors':errors,'per_species':per_species,
          'success_gate_necessary_only':checks,'immediate_stop_triggers':immediate,
          'decision':'NECESSARY_GATES_PASS_FURTHER_VALIDATION_REQUIRED' if all(checks.values()) else 'R1_NO_GO',
          'ablation_gates':'NOT_EXECUTED_NOT_B_MINUS_A','Setaria':'NOT_ACCESSED',
          'exact_coding_gene_F1_definition':'same as complete-primary exact chain F1 in frozen evaluator'})
    (OUT/'STATUS').write_text('COMPLETED_R1\n')
    print(json.dumps({'metrics':metrics,'gates':checks,'immediate_stop':immediate}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['infer','replay','decode','evaluate'])
    parser.add_argument('--arm',choices=['A','B'])
    parser.add_argument('--output-dir',type=Path,default=OUT)
    parser.add_argument('--score-source-run',type=Path,
                        help='Read saved scores from this run without copying or modifying them')
    args=parser.parse_args()
    OUT=args.output_dir
    SCORE_SOURCE=args.score_source_run
    if args.stage=='evaluate':
        if args.arm is None: parser.error('evaluate requires --arm')
        evaluate(args.arm)
    else:
        globals()[args.stage]()
