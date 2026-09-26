import itertools
import random
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION'))
import run_graph_core as G


def brute(sequence,exons,prefix,strand):
    genes=[]
    def visit(chain,score):
        e=chain[-1]
        if e.right=='stop':
            coding=''.join(sequence[x.a:x.b] for x in chain)
            if (len(coding)>=6 and len(coding)%3==0 and coding.startswith('ATG') and
                coding[-3:] in G.STOPS and not any(coding[i:i+3] in G.STOPS for i in range(3,len(coding)-3,3))):
                genes.append((tuple((x.a,x.b) for x in chain),score))
            return
        for f in exons:
            if f.left=='acceptor' and f.run>e.run and f.a>=e.b+4:
                visit(chain+[f],score+float(prefix[f.a]-prefix[e.b])+f.score)
    for e in exons:
        if e.left=='start': visit([e],e.score)
    # All complete paths and all compatible gene sets; independent small oracle.
    best=(0.,0,())
    def select(start,chosen,score,span):
        nonlocal best
        key=tuple(sorted(G.genomic_chain(g,len(sequence),strand) for g in chosen))
        value=(score,span,key)
        if (-value[0],value[1],value[2]) < (-best[0],best[1],best[2]): best=value
        for chain,s in genes:
            if chain[0][0]>=start:
                select(chain[-1][1],chosen+(chain,),score+s,span+chain[-1][1]-chain[0][0])
    select(0,(),0.,0)
    return best


def compare(sequence,exons,prefix=None):
    if prefix is None: prefix=np.zeros(len(sequence)+1)
    for strand in ('+','-'):
        expected=brute(sequence,exons,prefix,strand)
        actual=G.solve(sequence,exons,prefix,strand)
        key=tuple(sorted(G.genomic_chain(g,len(sequence),strand) for g in actual.genes))
        assert (actual.score,actual.span,key)==expected
    return G.solve(sequence,exons,prefix)


def test_cross_block_skip_long_intron():
    s='ATG'+'C'*20000+'AAATAA'
    exons=[G.Exon(0,0,3,'start','donor',4),G.Exon(4,20003,20009,'acceptor','stop',5)]
    p=compare(s,exons)
    assert p.genes==(((0,3),(20003,20009)),)
    assert G.models(p)[0]['phase']==[0,0]


def test_internal_join_stop():
    s='ATGT'+'C'*6+'AATAA'
    p=compare(s,[G.Exon(0,0,4,'start','donor',20),G.Exon(1,10,15,'acceptor','stop',20)])
    assert not p.genes


def test_short_initial_and_final_exons():
    # The actual ATG and final stop may straddle exon joins.
    s='ATGTCCCCAGGAAATAA'
    exons=[G.Exon(0,0,2,'start','donor',2),G.Exon(1,10,17,'acceptor','stop',2)]
    assert compare(s,exons).genes==(((0,2),(10,17)),)
    assert G.ORF(s).advance('',exons[0])=='INIT_AT'
    assert G.ORF('ATGAAATCCCCTAA').advance('T',G.Exon(1,12,14,'acceptor','stop',1))==''


def test_multiple_genes_compete_with_fusion():
    s='ATGAAATAA'+'C'*9+'ATGAAATAA'
    e=[G.Exon(0,0,9,'start','stop',7),G.Exon(2,18,27,'start','stop',7),
       G.Exon(0,0,6,'start','donor',4),G.Exon(2,21,27,'acceptor','stop',4)]
    p=compare(s,e)
    assert len(p.genes)==2 and p.score==14
    prefix=np.arange(len(s)+1,dtype=float)
    assert len(compare(s,e,prefix).genes)==1


def test_nonpositive_ties_run_order_and_coordinates():
    s='ATGAAATAA'+'C'*9+'ATGAAATAA'
    for score in (-1,0):
        assert not compare(s,[G.Exon(0,0,9,'start','stop',score)]).genes
    e=[G.Exon(0,0,9,'start','stop',2),G.Exon(1,0,27,'start','stop',2)]
    assert compare(s,e).span==9
    assert G.genomic_chain(((2,4),(10,15)),30,'-')==((15,20),(26,28))
    # Same-run and decreasing-run introns are disallowed even when geometric.
    for rid in (0,1):
        assert not compare(s,[G.Exon(1,0,3,'start','donor',5),G.Exon(rid,21,27,'acceptor','stop',5)]).genes


def test_suffix_automaton_exhaustive():
    # Independent concatenated-string check, including N and tiny terminal CDS.
    for state in G.STATES:
        for fragment in map(''.join,itertools.product('ACGTN',repeat=3)):
            sequence='CCC'+fragment
            e=G.Exon(1,3,6,'acceptor','donor',0)
            concat=state+fragment
            bad=any(concat[i:i+3] in G.STOPS for i in range(0,len(concat)-2,3))
            actual=G.ORF(sequence).advance(state,e)
            assert (actual is None)==bad
            if not bad:
                tail=concat[len(concat)//3*3:]
                assert actual==tail.replace('N','C')


def test_random_graphs_exhaustive():
    rng=random.Random(701)
    # Integer scores make tie-breaking exact, not tolerance-based.
    s='ATGAAATAA'+'CC'+'ATGCCCTAA'+'CC'+'ATGAAATAA'
    for trial in range(250):
        exons=[]
        for run,(a,b) in enumerate(((0,9),(11,20),(22,31))):
            for left,right,x,y in [('start','stop',a,b),('start','donor',a,a+3),
                                   ('acceptor','stop',a+3,b),('acceptor','donor',a+3,b-3)]:
                if rng.random()<0.8:
                    exons.append(G.Exon(run,x,y,left,right,rng.randrange(-3,4)))
        prefix=np.r_[0,np.cumsum([rng.randrange(-1,2) for _ in s])]
        compare(s,exons,prefix)


def test_candidate_threshold_equality():
    # Real candidate generation, including frozen inclusive threshold and motif.
    s='CCC'+'ATGAAATAA'+'CCC'
    r=np.zeros((len(s),3),dtype=np.float16)
    r[:,0]=4
    r[3:12,0]=0
    r[3:12,1]=4
    boundary=np.full((len(s),4),-10,dtype=np.float16)
    import redecode_error_decomposition as D
    boundary[3,D.m25.BOUNDARY_NAMES.index('start')]=0
    boundary[9,D.m25.BOUNDARY_NAMES.index('stop')]=0
    exons,prefix,counts=G.candidates(s,r,boundary)
    assert any((e.a,e.b,e.left,e.right)==(3,12,'start','stop') for e in exons)
    boundary[3,D.m25.BOUNDARY_NAMES.index('start')]=-0.001
    exons,_,_=G.candidates(s,r,boundary)
    assert not any(e.left=='start' for e in exons)
