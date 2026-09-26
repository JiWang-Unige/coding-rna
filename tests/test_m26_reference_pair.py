import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import reference_pair as R


def test_araport_mapping_keeps_coordinates_and_does_not_invent_coding_genes():
    coding = 'Chr3\tAraport11\tgene\t10\t20\t.\t-\t.\tID=AT3G00010;locus_type=protein_coding\n'
    other = coding.replace('protein_coding', 'lncRNA')
    assert R.normalize_araport(coding).startswith(R.SEQID+'\tAraport11\tgene\t10\t20')
    assert 'gene_biotype=protein_coding' in R.normalize_araport(coding)
    assert 'gene_biotype=lncRNA' in R.normalize_araport(other)
    assert R.normalize_araport(coding.replace('Chr3','Chr4')) is None


def test_longest_CDS_primary_and_transcript_ID_tie_break():
    txs = [{'id':'b','CDS':[(1,10,'0')]}, {'id':'a','CDS':[(1,10,'0')]}, {'id':'c','CDS':[(1,7,'0')]}]
    assert R.pick_primary({'g':txs})[0]['id'] == 'a'
