import sys
from pathlib import Path
import numpy as np
import pytest

BASE=Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY"
sys.path.insert(0,str(BASE))
import confirm_run as C
import wt_tiberius as T
from test_m27_paired_response import example


def test_ids_ignored_but_different_bypass_structure_does_not_pass():
    row,dna,tx=example()
    wt={C.key(tx["WT"])}
    old=[tx["bypass"]]
    renamed=[{**tx["bypass"],"id":"different"}]
    assert C.compare_target(row,old,renamed,dna,wt,{})["equal"]
    other={**tx["bypass"],"CDS":[(0,9,0),(20,23,0),(26,29,0),(40,49,0)]}
    result=C.compare_target(row,old,[other],dna,wt,{})
    assert result["old"]["category"]==result["new"]["category"]=="BYPASS"
    assert not result["equal"]


def test_multiplicity_phase_and_competition_are_not_ignored():
    row,dna,tx=example()
    old=[tx["bypass"]]
    wt={C.key(tx["WT"])}
    for new in (old+[{**tx["bypass"],"id":"duplicate"}],old+[tx["prefix"]],
                [{**tx["bypass"],"CDS":[(0,9,0),(40,49,1)]}]):
        assert not C.compare_target(row,old,new,dna,wt,{})["equal"]


def test_non_target_differences_are_saved_without_failing_target_identity():
    row,dna,tx=example()
    old=[tx["bypass"]]
    new=old+[{"id":"remote","seqid":"chr22","strand":"+","CDS":[(54,57,0)]}]
    assert C.compare_target(row,old,new,dna,{C.key(tx["WT"])},{})["equal"]
    diff=C.chain_difference(old,new)
    assert not diff["equal"] and diff["lost"]==[] and diff["added"][0]["count"]==1


def test_filter_target_change_detected_even_when_final_chain_same():
    row,dna,tx=example()
    bypass=[C.key(tx["bypass"])]
    original=[C.key(tx["WT"])]
    names=["check_min_coding_length","check_inframe_stop_codons"]
    old={"filters":[{"name":n,"before":bypass,"after":bypass} for n in names]}
    new={"filters":[{"name":names[0],"before":original,"after":bypass},old["filters"][1]]}
    result=C.compare_filters(row,old,new,dna,{C.key(tx["WT"])},{})
    assert not result[0]["before"]["target"]["equal"]
    assert result[0]["after"]["target"]["equal"]


def test_record_keeps_native_object_and_budget_includes_file_headers(tmp_path):
    x=np.zeros((1,3,5),dtype=np.float32)
    y=np.ones((1,3,15),dtype=np.float32)
    small=T.NeuralCache(tmp_path/"small","record",limit=x.nbytes+y.nbytes)
    with pytest.raises(RuntimeError,match="limit"):
        small.call(x,{},lambda:y)
    assert not list((tmp_path/"small").glob("*.npy"))
    cache=T.NeuralCache(tmp_path/"cache","record",limit=8192)
    assert cache.call(x,{"stage":"initial"},lambda:y) is y
    cache.finish()
    assert 0<cache.bytes<8192
