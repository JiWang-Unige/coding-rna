import importlib.util
from pathlib import Path

path = Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY/prepare_reference.py"
spec = importlib.util.spec_from_file_location("m27_prepare",path)
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def fixture(strand="+"):
    seq = "ATG"+"GAA"*39+"TAT"+"AAA"+"GAA"*48
    genome = list("N"*500)
    intervals = [(10,100),(160,250),(310,400)]
    for p,b in zip(M.ordered_positions(intervals,"+"),seq):
        genome[p]=b
    genome[400:403]=list("TGA")
    tx={"attrs":{},"strand":"+","CDS":[(s,e,0) for s,e in intervals],
        "stop_codon":[(400,403)],"start_codon":[(10,13)],"selenocysteine":False}
    if strand == "-":
        genome=list("".join(genome).translate(M.COMP)[::-1])
        tx["strand"]="-"
        tx["CDS"]=[(500-e,500-s,ph) for s,e,ph in tx["CDS"]]
        tx["stop_codon"]=[(97,100)]
        tx["start_codon"]=[(487,490)]
    return tx,"".join(genome)


def test_plus_and_minus_same_base_pair_and_stop_union():
    for strand in ("+","-"):
        tx,genome=fixture(strand)
        obj,reason=M.assess_transcript(tx,genome)
        assert reason=="assessable"
        assert len(obj["seq"])==273 and obj["seq"][-3:]=="TGA"
        site=M.paired_site(tx,obj,genome)
        assert site is not None
        q=site["position_1based"]-1
        assert genome[q]==site["reference_base"]
        for field,codon in (("synonymous_alt","TAC"),("stop_alt","TAA")):
            mutated=genome[:q]+site[field]+genome[q+1:]
            new=M.bases(mutated,obj["positions"],strand)
            assert new[120:123]==codon
            assert sum(a!=b for a,b in zip(mutated,genome))==1


def test_internal_stop_and_phase_are_not_assessable():
    tx,genome=fixture()
    q=M.ordered_positions([(s,e) for s,e,_ in tx["CDS"]],"+")[122]
    altered=genome[:q]+"A"+genome[q+1:]
    assert M.assess_transcript(tx,altered)[1]=="internal_stop"
    tx["CDS"][1]=(160,250,1)
    assert M.assess_transcript(tx,genome)[1]=="CDS_phase"


def test_stop_feature_must_not_overlap_CDS():
    tx,genome=fixture()
    tx["stop_codon"]=[(397,400)]
    assert M.assess_transcript(tx,genome)[1]=="CDS_stop_geometry"


def test_new_canonical_motif_excluded():
    tx,genome=fixture()
    obj,_=M.assess_transcript(tx,genome)
    # TAT at 190..192; next base G creates AG if Tyr third base changes to A.
    q=obj["positions"][122]
    genome=genome[:q+1]+"G"+genome[q+2:]
    assert M.paired_site(tx,obj,genome) is None


def test_quantile_selection_reference_only_and_no_duplicates():
    rows=[{"gene_id":str(i),"position_1based":i} for i in range(101)]
    chosen=M.choose_loci(rows)
    assert len(chosen)==24 and len({r["gene_id"] for r in chosen})==24
    assert M.choose_loci(rows[:3])==rows[:3]
    assert M.choose_loci([])==[]
