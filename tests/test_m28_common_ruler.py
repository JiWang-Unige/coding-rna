import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import common_ruler as E


def test_assessable_primary_does_not_inherit_parent_partial_or_training_geometry():
    dna = 'ATGAAATAA' + 'N'*20 + 'ATGAAAA'
    good = {'id':'good', 'gene_id':'g', 'seqid':'q', 'strand':'+', 'CDS':[(0,9,'0')],
            'partial':True, 'exon':[]}
    bad = dict(good, id='longest_bad', CDS=[(0,len(dna),'0')])
    refs, all_tx = E.select_reference({'transcripts':{'good':good,'bad':bad}}, {}, {'q':dna})
    assert refs == all_tx == [good]
    assert not E.select_reference({'transcripts':{'good':good}}, {'good':[{'partial':'true'}]}, {'q':dna})[0]


def test_other_isoform_and_all_gene_background_are_separate_from_primary_truth():
    def tx(i,a,b):
        return {'id':i,'gene_id':i,'seqid':'q','strand':'+','CDS':[(a,b,'0')]}
    primary, iso, outside = tx('p',0,9), tx('i',0,12), tx('b',30,39)
    metrics = E.score([iso,outside,outside], [primary], [primary,iso], 50,
                      {'genes':{'coding':{'start':0,'end':15},'noncoding':{'start':20,'end':25}}})
    assert metrics['exact_chain']['tp'] == 0
    assert metrics['exact_other_assessable_isoform_not_primary_TP'] == 1
    assert metrics['duplicate_chain_records'] == 1
    assert metrics['background']['bp'] == 30
    assert metrics['background']['wholly_background_unique_chains'] == 1
    assert metrics['background']['predicted_span_bp'] == 9
