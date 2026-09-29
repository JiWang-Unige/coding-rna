import sys
from pathlib import Path
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/'scripts')]
from src.m28.export import as_transcript,write_gff3
import eval_structure_diagnostic as E

def test_primary_gff_roundtrip_halfopen_coordinates_and_both_strand_phases(tmp_path):
    chain=[(100,104),(120,128)]
    plus=as_transcript(chain,'q','+','p',1.25)
    minus=as_transcript(chain,'q','-','m',2.)
    assert plus['CDS']==[(100,104,'0'),(120,128,'2')]
    assert minus['CDS']==[(100,104,'1'),(120,128,'0')]
    path=tmp_path/'pred.gff3'
    write_gff3(path,[plus,minus],{'q':200},'M28')
    parsed=E.parse_annotation(path,{'q':200},protein_coding_only=True)
    assert len(parsed['genes'])==2
    for t in (plus,minus):
        actual=parsed['transcripts'][t['id']]
        assert actual['CDS']==t['CDS']
        assert actual['strand']==t['strand']
    text=path.read_text()
    assert '\tCDS\t101\t104\t' in text
    assert 'UTR' not in text

def test_empty_gff_remains_parseable(tmp_path):
    path=tmp_path/'empty.gff3'
    write_gff3(path,[],{'q':200},'M28')
    assert not E.parse_annotation(path,{'q':200})['transcripts']
