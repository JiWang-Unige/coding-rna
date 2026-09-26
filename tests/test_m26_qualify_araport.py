import io
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import qualify_araport as Q


def test_type_comes_from_same_release_table_not_free_text_or_CDS_presence():
    attrs = Q.structural_attrs('ID=AT3G01560;Dbxref=x,locus_type equals protein_coding;Note=unescaped;words')
    assert Q.decide_type('AT3G01560', attrs, {}) == (None, 'unresolved')
    assert Q.decide_type('AT3G01560', attrs, {'AT3G01560': {'protein_coding'}}) == ('protein_coding', 'same_release_functional_table')
    assert Q.decide_type('AT3G01560', attrs, {'AT3G01560': {'small_nuclear_rna'}})[0] == 'small_nuclear_rna'


def test_conflicting_gene_types_block_qualification():
    assert Q.decide_type('AT3G01560', {'locus_type': 'protein_coding'}, {'AT3G01560': {'pseudogene'}})[0] is None
    assert Q.decide_type('AT3G01560', {}, {'AT3G01560': {'pseudogene', 'protein_coding'}})[0] is None


def test_table_preserves_isoform_disagreement_and_structural_attributes():
    mapping, _, _, _ = Q.type_table(iter(io.StringIO('name\tgene_model_type\tdescription\nAT3G01560.1\tprotein_coding\ttext\nAT3G01560.2\tpseudogene\ttext\n')))
    assert len(mapping['AT3G01560']) == 2
    attrs = Q.structural_attrs('ID=x;Parent=g;partial=true;start_range=.,42;Note=a;b;transl_except=example')
    assert attrs == {'ID':'x', 'Parent':'g', 'partial':'true', 'start_range':'.,42', 'transl_except':'example'}


def test_real_uORF_records_are_retained_without_overwriting_host_gene_type():
    text = 'name\tgene_model_type\tdescription\nAT3G25570.1\tprotein_coding\tgene\nAT3G25570.uORF1-2\tuORF\tDerives_from AT3G25570;(source:Araport11)\n'
    mapping, _, subfeatures, unknown = Q.type_table(iter(io.StringIO(text)))
    assert mapping['AT3G25570'] == {'protein_coding'}
    assert subfeatures[0]['model'] == 'AT3G25570.uORF1-2'
    assert unknown == []
    _, _, _, unknown = Q.type_table(iter(io.StringIO(text.replace('Derives_from AT3G25570', 'Derives_from AT3G59050'))))
    assert unknown[0]['reason'] == 'uORF_host_not_confirmed'
