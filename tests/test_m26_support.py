import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import support_partition as P


def test_any_ordered_witness_not_only_historical_witness():
    assert P.witness([[8,2],[4,9],[6,10]]) == [2,4,6]
    assert P.witness([[2],[1,2]]) is None
    assert P.witness([[3],[7,9],[8]]) == [3,7,8]


def test_actual_membership_requires_grammar_and_filtered_fragments():
    orf = P.G.ORF('ATGAAATAA')
    assert P.membership([(0,9)],{(0,9,'start','stop'):[2]},orf)
    assert not P.membership([(0,9)],{},orf)
    assert not P.membership([(0,9)],{(0,9,'start','stop'):[2]},P.G.ORF('ATGAAACCC'))


def test_simultaneous_bound_accounts_for_same_strand_conflicts():
    def tx(a,b,strand):
        return {'seqid':'q','strand':strand,'CDS':[(a,b,'0')]}
    assert P.simultaneous_bound([tx(0,10,'+'),tx(2,5,'+'),tx(5,8,'+'),tx(0,10,'-')]) == 3
