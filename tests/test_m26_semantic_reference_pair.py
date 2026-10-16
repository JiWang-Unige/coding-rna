import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import semantic_reference_pair as R


def test_primary_selection_is_not_structural_absence():
    a, b = ('chr','+',((1,10),)), ('chr','+',((1,20),))
    assert R.locus_primary_class(a, b, {a,b}, {a,b}) == 'primary_selection_only'
    assert R.locus_primary_class(a, b, {a}, {b}) == 'primary_structure_disagreement'
    assert R.locus_primary_class(a, a, {a,b}, {a}) == 'same_primary'


def test_one_sided_assessability_is_not_reported_as_competing_primary_structures():
    a = ('chr','+',((1,10),))
    assert R.locus_primary_class(None, a, set(), {a}) == 'Araport_only_CDS_assessable'
    assert R.locus_primary_class(a, None, {a}, set()) == 'RefSeq_only_CDS_assessable'
    assert R.locus_primary_class(None, None, set(), set()) == 'neither_CDS_assessable'


def test_only_observed_explicit_miRNA_vocabulary_difference_is_noncoding():
    attrs = {'locus_type':'mirna'}
    assert R.classify_type('g', 'gene', attrs, {'g':{'miRNA_primary_transcript'}}) == ('mirna','explicit_noncoding_vocabulary_disagreement')
    assert R.classify_type('g', 'gene', attrs, {'g':{'protein_coding'}})[0] is None


def test_TE_or_unknown_cannot_be_promoted_to_protein_coding():
    assert R.classify_type('g', 'transposable_element_gene', {'locus_type':'protein_coding'}, {'g':{'protein_coding'}})[0] is None
    assert R.classify_type('g', 'gene', {}, {})[0] is None
    assert R.classify_type('g', 'transposable_element_gene', {'locus_type':'transposable_element_gene'}, {'g':{'transposable_element_gene'}})[0] == 'transposable_element_gene'


def test_CDS_assessability_difference_is_not_structure_absence():
    a, b = ('chr','+',((1,10),)), ('chr','+',((1,20),))
    assert R.locus_primary_class(a, b, {a}, {b}, {a,b}, {a,b}) == 'CDS_assessability_difference'


def test_phase_chain_is_checked_in_transcription_order():
    plus = {'strand':'+','CDS':[(1,5,'0'),(10,15,'2')]}
    minus = {'strand':'-','CDS':[(1,6,'2'),(10,14,'0')]}
    assert R.phase_chain_consistent(plus)
    assert R.phase_chain_consistent(minus)
    assert not R.phase_chain_consistent({'strand':'+','CDS':[(1,5,'0'),(10,15,'1')]})


def test_explicit_pseudogene_vs_pre_tRNA_stays_noncoding_without_claiming_equivalence():
    assert R.classify_type('AT3G20365', 'pseudogene', {'locus_type':'pseudogene'}, {'AT3G20365':{'pre_trna'}}) == ('pseudogene','explicit_noncoding_pseudogene_tRNA_disagreement')
    assert R.classify_type('AT3G20365', 'pseudogene', {'locus_type':'pseudogene'}, {'AT3G20365':{'protein_coding'}})[0] is None
