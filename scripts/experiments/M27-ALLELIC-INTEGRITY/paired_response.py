#!/usr/bin/env python3
"""Frozen, conservative M27 target-response classification (not biological truth)."""
import argparse
import json
from collections import Counter
from pathlib import Path
from prepare_reference import attributes, bases, merged, ordered_positions, STOP
from wt_compare import key, reference_key

ROOT=Path(__file__).resolve().parents[3]


def subtract(intervals, masks):
    result=[]
    for start,end in merged(intervals):
        cursor=start
        for a,b in merged(masks):
            if b<=cursor:
                continue
            if a>=end:
                break
            if cursor<a:
                result.append((cursor,min(a,end)))
            cursor=max(cursor,b)
            if cursor>=end:
                break
        if cursor<end:
            result.append((cursor,end))
    return result


def intersects(a,b):
    return max(a[0],b[0])<min(a[1],b[1])


def reference_coding(path):
    genes={}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f=line.rstrip().split("\t")
            if f[0]!="chr22" or f[2] not in {"CDS","stop_codon"}:
                continue
            gid=attributes(f[8])["gene_id"][0]
            item=genes.setdefault(gid,{"strand":f[6],"intervals":[]})
            item["intervals"].append((int(f[3])-1,int(f[4])))
    return {g:{"strand":r["strand"],"intervals":merged(r["intervals"])} for g,r in genes.items()}


def foreign_regions(row,genes):
    own=genes[row["gene_id"]]["intervals"]
    return {g:subtract(r["intervals"],own) for g,r in genes.items()
            if g!=row["gene_id"] and r["strand"]==row["strand"]}


def blocks(tx):
    return sorted((s,e,int(p)) for s,e,p in tx["CDS"])


def ordered(tx):
    return sorted(blocks(tx),reverse=tx["strand"]=="-")


def positions(tx):
    return ordered_positions([(s,e) for s,e,_ in tx["CDS"]],tx["strand"])


def introns(tx):
    parts=blocks(tx)
    return [(a[1],b[0]) for a,b in zip(parts,parts[1:])]


def coding_status(tx,genome):
    cumulative=0
    phase_ok=True
    for s,e,p in ordered(tx):
        phase_ok=phase_ok and p==(3-cumulative%3)%3
        cumulative+=e-s
    sequence=bases(genome,positions(tx),tx["strand"])
    canonical=(len(sequence)>=6 and len(sequence)%3==0 and not set(sequence)-set("ACGT")
               and sequence[:3]=="ATG" and sequence[-3:] in STOP)
    internal=sum(sequence[i:i+3] in STOP for i in range(0,max(0,len(sequence)-3),3))
    return {"phase_continuous":phase_ok,"canonical_start_terminal_stop":canonical,
            "internal_stop_count":internal,"valid_ORF":bool(phase_ok and canonical and internal==0)}


def classify(row,predictions,genome,wt_keys,foreign):
    refkey=reference_key(row)
    ref={"seqid":"chr22","strand":row["strand"],"CDS":list(refkey[2])}
    refpos=positions(ref)
    codon_end=3*row["codon_index_1based"]
    codon=refpos[codon_end-3:codon_end]
    first,last=ordered(ref)[0],ordered(ref)[-1]
    span=(refkey[2][0][0],refkey[2][-1][1])
    related=[t for t in predictions if t["seqid"]=="chr22" and t["strand"]==row["strand"] and t["CDS"]
             and intersects((min(s for s,_,_ in t["CDS"]),max(e for _,e,_ in t["CDS"])),span)]
    result={"related_prediction_ids":[t["id"] for t in related],
            "WT_first_CDS_anchor":list(first),"WT_last_CDS_anchor":list(last),
            "edited_codon_positions_0based":codon}
    def done(category,bounds,assigned=(),**extra):
        return {**result,"category":category,"bypass_interval":list(bounds),
                "assigned_prediction_ids":[t["id"] for t in assigned],
                "assigned_keys":[key(t) for t in assigned],**extra}
    if not related:
        return done("OUTPUT_ABSENT",(0,0))
    dual=[t for t in related if ordered(t)[0]==first and ordered(t)[-1]==last]
    if len(dual)>1:
        return done("AMBIGUOUS_CANDIDATES",(0,1),reason="multiple_double_anchor_chains")
    if len(dual)==1:
        tx=dual[0]
        competing=[t["id"] for t in related if t is not tx and (first in blocks(t) or last in blocks(t))]
        if competing:
            return done("UNASSIGNABLE",(0,1),reason="competing_terminal_anchor",competitors=competing)
        added=subtract([(s,e) for s,e,_ in tx["CDS"]],row["CDS_chain_0based_halfopen"])
        foreign_hits=[g for g,regions in foreign.items() if any(intersects(a,b) for a in added for b in regions)]
        if foreign_hits:
            return done("FUSION_SUSPECT",(0,1),foreign_coding_genes=foreign_hits)
        status=coding_status(tx,genome)
        if key(tx)==refkey:
            return done("WT_CHAIN_RETAINED",(0,0),(tx,),coding_status=status)
        excludes=any(all(a<=q<b for q in codon) for a,b in introns(tx))
        changed=introns(tx)!=introns(ref)
        if excludes and changed and status["valid_ORF"]:
            return done("BYPASS",(1,1),(tx,),coding_status=status)
        return done("OTHER_IDENTIFIED_STRUCTURE",(0,0),(tx,),coding_status=status,
                    codon_entirely_intronic=excludes,splice_connectivity_changed=changed)
    prefix=[]
    restart=[]
    upstream=[]
    for tx in related:
        pos=positions(tx)
        status=coding_status(tx,genome)
        is_background=key(tx) in wt_keys
        if (pos==refpos[:codon_end] and status["valid_ORF"]
                and bases(genome,codon,row["strand"]) in STOP):
            prefix.append(tx)
        if pos and pos[0] in refpos:
            j=refpos.index(pos[0])
            if j>=codon_end and pos==refpos[j:] and status["valid_ORF"] and not is_background:
                restart.append(tx)
        if ((first[1]-first[0])<=len(pos)<=codon_end
                and pos==refpos[:len(pos)] and status["phase_continuous"] and not is_background):
            upstream.append(tx)
    if len(related)==1 and len(prefix)==1:
        return done("PTC_TERMINATION",(0,0),prefix)
    if len(related)==1 and len(restart)==1:
        return done("DOWNSTREAM_RESTART",(0,0),restart)
    if (len(related)==2 and len(upstream)==1 and len(restart)==1
            and upstream[0] is not restart[0]):
        return done("OUTPUT_SPLIT",(0,0),(upstream[0],restart[0]))
    return done("UNASSIGNABLE",(0,1),reason="overlapping_output_without_unique_frozen_identity")


def snapshot_txs(records):
    return [{"id":f"trace_{i}","seqid":seq,"strand":strand,"CDS":[tuple(b) for b in parts]}
            for i,(seq,strand,parts) in enumerate(records)]


def filter_labels(row,trace,genome,wt_keys,foreign):
    result=[]
    for event in trace["filters"]:
        before=snapshot_txs(event["before"])
        after=snapshot_txs(event["after"])
        interpretation=classify(row,before,genome,wt_keys,foreign)
        assigned=interpretation["assigned_keys"]
        after_keys=Counter(map(key,after))
        removed=(len(assigned)==1 and after_keys[assigned[0]]==0)
        result.append({"filter":event["name"],"before_category":interpretation["category"],
                       "after_category":classify(row,after,genome,wt_keys,foreign)["category"],
                       "removed_unique_identified_chain":removed,
                       "removed_key":assigned[0] if removed else None})
    return result


def paired_bounds(pairs):
    n=len(pairs)
    low=sum(p["PTC"]["bypass_interval"][0]-p["syn"]["bypass_interval"][1] for p in pairs)
    high=sum(p["PTC"]["bypass_interval"][1]-p["syn"]["bypass_interval"][0] for p in pairs)
    table=Counter()
    for p in pairs:
        a,b=p["PTC"]["bypass_interval"],p["syn"]["bypass_interval"]
        table[f"{a[0]}{b[0]}" if a[0]==a[1] and b[0]==b[1] else "unknown_pair"]+=1
    return {"fixed_denominator":n,"difference_lower":low/n,"difference_upper":high/n,
            "bound_type":"identification_bounds_not_confidence_interval",
            "paired_table_PTC_syn":{k:table[k] for k in ("00","01","10","11","unknown_pair")},
            "all_pairs_identified":table["unknown_pair"]==0,"positive_lower_bound":low>0}


def summarize(out):
    manifest=json.loads((out/"input_manifest.json").read_text())
    rows=[]
    completed=0
    for row in manifest["registered_loci"]:
        record={"ordinal":row["ordinal"],"gene_id":row["gene_id"],"gene_name":row["gene_name"],
                "position_1based":row["position_1based"]}
        for method,state in row["method_status"].items():
            if state!="paired_native_planned":
                record[method]={"WT_status":"WT_NOT_EXACT","conditional_response":"NOT_RUN"}
                continue
            record[method]={"WT_status":"EXACT","conditional_response":"PLANNED"}
            for allele in ("syn","PTC"):
                base=out/method/f'{row["ordinal"]:02d}_{allele}'
                if ((base/"response.json").is_file() and (base/"execution.json").is_file()
                        and json.loads((base/"execution.json").read_text())["status"]=="COMPLETED"):
                    record[method][allele]=json.loads((base/"response.json").read_text())
                    completed+=1
                else:
                    record[method][allele]={"category":"TECHNICALLY_UNCOMPLETED","bypass_interval":[0,1]}
        rows.append(record)
    groups={}
    common=[r for r in rows if all(r[m]["WT_status"]=="EXACT" for m in ("ANNEVO","Tiberius"))]
    for m in ("ANNEVO","Tiberius"):
        eligible=[r for r in rows if r[m]["WT_status"]=="EXACT"]
        candidates=[r for r in eligible if r[m]["PTC"]["bypass_interval"]==[1,1] and r[m]["syn"]["bypass_interval"]==[0,0]]
        groups[m]={**paired_bounds([r[m] for r in eligible]),
                   "common_WT":paired_bounds([r[m] for r in common]),
                   "first_preregistered_mechanism_candidate":candidates[0]["gene_id"] if candidates else None}
    result={"stage":"native_paired_response_not_component_attribution","planned_method_cases":44,
            "completed_classified_cases":completed,"technically_complete":completed==44,
            "registered_loci_count":len(rows),"common_WT_count":len(common),"methods":groups,"loci":rows,
            "training_runs":0,"Setaria_accessed":False}
    with (out/"paired_result.json").open("w") as handle:
        json.dump(result,handle,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!="loci"},indent=2))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    summarize(args.output_dir)
