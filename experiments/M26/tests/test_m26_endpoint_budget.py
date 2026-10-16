import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import endpoint_budget as E


def test_rank_budget_threshold_tie_and_window_sham():
    pos=np.array([1,3,7,9,12,15])
    scores=np.array([0.,3.,3.,-4.,1.,2.])
    assert E.select(pos,scores,3,0).tolist()==[3,7,15]
    shuffled=E.permute_windows(pos,scores,np.random.default_rng(99),8)
    assert sorted(shuffled[:3])==sorted(scores[:3])
    assert sorted(shuffled[3:])==sorted(scores[3:])
    assert np.array_equal(shuffled,E.permute_windows(pos,scores,np.random.default_rng(99),8))
    assert len(E.select(pos,shuffled,3,0))==3


def test_endpoint_conversion_and_membership_without_run_constraint():
    s='ATGGTCCCCAGAAATAA'
    chain=[(0,3),(11,17)]
    exons=[E.G.Exon(8,0,3,'start','donor',0),E.G.Exon(2,11,17,'acceptor','stop',0)]
    points=E.endpoints(exons)
    assert points['acceptor'].tolist()==[9] and points['stop'].tolist()==[14]
    assert all(E.exists(points[e],p) for e,p in E.events_for(chain))
    assert E.grammar(chain,E.G.ORF(s))
    assert not E.P.membership(chain,{(x.a,x.b,x.left,x.right):[x.run] for x in exons},E.G.ORF(s))
    # A reverse-strand chain must use transcription-oriented coordinates.
    assert E.G.genomic_chain(tuple(chain),len(s),'-')==((0,6),(14,17))


def test_motifs_and_ORF_failures():
    s='ATGGTCCCCAGAAATAA'
    assert E.motifs(s,'start').tolist()==[0]
    assert E.motifs(s,'donor').tolist()==[3]
    assert E.motifs(s,'acceptor').tolist()==[9]
    assert E.motifs(s,'stop').tolist()==[14]
    assert not E.grammar([(0,3),(6,9)],E.G.ORF('ATGCCCTAA'))
    assert not E.grammar([(0,12)],E.G.ORF('ATGTAACCCTAA'))


def test_pruned_fragment_enumeration_equals_exhaustive_with_short_CDS():
    # All motif endpoints and all endpoint-pairs; upper_end is merely an accelerator.
    for s in ('ATGGTCCCCAGAAATAA','ATGATGTAGTAAAGTGACCTAA','ATGTCCCCAGGAAATAA'):
        orf=E.G.ORF(s)
        points={e:E.motifs(s,e) for e in E.EVENTS}
        expected=set()
        for l in ('start','acceptor'):
            for r in ('donor','stop'):
                for p in points[l]:
                    for q in points[r]:
                        a,b=int(p)+(2 if l=='acceptor' else 0),int(q)+(3 if r=='stop' else 0)
                        x=E.G.Exon(0,a,b,l,r,0.)
                        if a<b and E.viable(x,orf):
                            expected.add((a,b,l,r))
        actual={(x.a,x.b,x.left,x.right) for x in E.flat_fragments(points,orf)}
        assert actual==expected
        counts=E.opportunity_counts(points,orf)
        assert counts['local_ORF_admissible_fragments']==len(expected)
        joins=sum(x[1]+4<=y[0] for x in expected for y in expected if x[3]=='donor' and y[2]=='acceptor')
        assert counts['geometric_intron_vertex_pair_upper_bound']==joins


def test_original_connection_count_enforces_run_order():
    s='ATGGTCCCCAGAAATAA'
    first=E.G.Exon(1,0,3,'start','donor',0)
    for rid,expected in ((0,0),(1,0),(2,1)):
        exons=[first,E.G.Exon(rid,11,17,'acceptor','stop',0)]
        counts=E.opportunity_counts(E.endpoints(exons),E.G.ORF(s),exons)
        assert counts['local_ORF_admissible_fragments']==2
        assert counts['geometric_intron_vertex_pair_upper_bound']==expected
