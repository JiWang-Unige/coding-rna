import importlib.util
import sys
from pathlib import Path
import pytest
BASE=Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY"
sys.path.insert(0,str(BASE))
import paired_response as M
from prepare_reference import COMP


def example(strand="+",allele="PTC"):
    text=list("N"*60)
    for start,seq in [(0,"ATGAAACCC"),(20,"AAATACATG"),(40,"CCCAAATAA")]:
        text[start:start+len(seq)]=seq
    wt="".join(text)
    row={"gene_id":"target","strand":"+","CDS_chain_0based_halfopen":[[0,9],[20,29],[40,49]],
         "codon_index_1based":5}
    point=25
    dna=wt[:point]+("A" if allele=="PTC" else "T")+wt[point+1:]
    def tx(name,parts):
        if strand=="-":
            parts=sorted((60-e,60-s,p) for s,e,p in parts)
        return {"id":name,"seqid":"chr22","strand":strand,"CDS":parts}
    if strand=="-":
        dna=dna.translate(COMP)[::-1]
        row["CDS_chain_0based_halfopen"]=[[60-e,60-s] for s,e in row["CDS_chain_0based_halfopen"][::-1]]
        row["strand"]="-"
    pred={
        "WT":tx("WT",[(0,9,0),(20,29,0),(40,49,0)]),
        "bypass":tx("bypass",[(0,9,0),(40,49,0)]),
        "prefix":tx("prefix",[(0,9,0),(20,26,0)]),
        "restart":tx("restart",[(26,29,0),(40,49,0)]),
        "readthrough":tx("readthrough",[(0,9,0),(40,55,0)]),
        "overlap":tx("overlap",[(10,19,0)]),
        "fusion":tx("fusion",[(0,9,0),(30,33,0),(40,49,0)]),
    }
    return row,dna,pred


@pytest.mark.parametrize("strand",["+","-"])
def test_strict_bypass_and_ptc_termination_on_both_strands(strand):
    row,dna,t=example(strand)
    keys={M.key(t["WT"])}
    assert M.classify(row,[t["bypass"]],dna,keys,{})["category"]=="BYPASS"
    assert M.classify(row,[t["prefix"]],dna,keys,{})["category"]=="PTC_TERMINATION"
    assert M.classify(row,[t["restart"]],dna,keys,{})["category"]=="DOWNSTREAM_RESTART"
    assert M.classify(row,[t["prefix"],t["restart"]],dna,keys,{})["category"]=="OUTPUT_SPLIT"


def test_background_restart_competing_anchor_and_readthrough_remain_unknown():
    row,dna,t=example()
    keys={M.key(t["WT"]),M.key(t["restart"])}
    for predictions in ([t["restart"]],[t["bypass"],t["prefix"]],[t["bypass"],{**t["bypass"],"id":"duplicate"}],[t["readthrough"]]):
        assert M.classify(row,predictions,dna,keys,{})["bypass_interval"]==[0,1]


def test_absence_is_not_missing_anchor_and_retained_internal_stop_is_not_repaired():
    row,dna,t=example()
    keys={M.key(t["WT"])}
    assert M.classify(row,[],dna,keys,{})["category"]=="OUTPUT_ABSENT"
    assert M.classify(row,[t["overlap"]],dna,keys,{})["category"]=="UNASSIGNABLE"
    r=M.classify(row,[t["WT"]],dna,keys,{})
    assert r["category"]=="WT_CHAIN_RETAINED" and r["coding_status"]["internal_stop_count"]==1


def test_foreign_coding_attachment_and_invalid_phase_do_not_become_bypass():
    row,dna,t=example()
    dna=dna[:30]+"AAA"+dna[33:]
    keys={M.key(t["WT"])}
    assert M.classify(row,[t["fusion"]],dna,keys,{"foreign":[(30,33)]})["category"]=="FUSION_SUSPECT"
    bad={**t["fusion"],"CDS":[(0,9,0),(30,33,1),(40,49,0)]}
    r=M.classify(row,[bad],dna,keys,{})
    assert r["category"]=="OTHER_IDENTIFIED_STRUCTURE" and not r["coding_status"]["phase_continuous"]


def test_filter_removal_is_observed_separately_from_final_endpoint():
    row,dna,t=example()
    def snap(tx):
        return [tx["seqid"],tx["strand"],tx["CDS"]]
    trace={"filters":[{"name":"check_min_coding_length","before":[snap(t["prefix"])],"after":[]},
                      {"name":"check_inframe_stop_codons","before":[],"after":[]}]}
    labels=M.filter_labels(row,trace,dna,{M.key(t["WT"])},{})
    assert labels[0]["removed_unique_identified_chain"]
    assert not labels[1]["removed_unique_identified_chain"]
    assert M.classify(row,[],dna,{M.key(t["WT"])},{})["bypass_interval"]==[0,0]


def test_fixed_denominator_identification_bounds_and_opposing_pairs():
    def state(a,b=None):
        return {"bypass_interval":[a,a if b is None else b]}
    pairs=[{"PTC":state(1),"syn":state(0)},{"PTC":state(0),"syn":state(1)},
           {"PTC":state(0,1),"syn":state(0)}]
    r=M.paired_bounds(pairs)
    assert r["fixed_denominator"]==3
    assert r["difference_lower"]==0 and r["difference_upper"]==pytest.approx(1/3)
    assert r["paired_table_PTC_syn"]=={"00":0,"01":1,"10":1,"11":0,"unknown_pair":1}


def test_exclusive_foreign_regions_do_not_count_target_shared_annotation():
    row,_,_=example()
    genes={"target":{"strand":"+","intervals":[(0,9),(20,29),(40,49)]},
           "other":{"strand":"+","intervals":[(20,29),(30,33)]},
           "antisense":{"strand":"-","intervals":[(30,33)]}}
    assert M.foreign_regions(row,genes)=={"other":[(30,33)]}
