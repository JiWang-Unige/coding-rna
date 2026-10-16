import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION'))
import run_graph_run as R
import run_graph_core as G


def test_split_codon_export_and_independent_audit(tmp_path,monkeypatch):
    monkeypatch.setattr(R,'OUT',tmp_path)
    sequence='ATGTCCCCAGGAAATAA'
    chain=((0,2),(10,17))
    model=G.models(G.Path(score=1,span=17,genes=(chain,)))[0]
    predictions={('synthetic','plus'):[R.D.m25._map_model_from_orientation(model,17,'+')],
                 ('synthetic','minus'):[R.D.m25._map_model_from_orientation(model,17,'-')]}
    R.export('B',predictions)
    seqs={'plus':sequence,'minus':R.D.m25.reverse_complement(sequence)}
    audit=R.audited_transcripts(tmp_path/'B_predictions.gff3',seqs,2,split_codons=True)
    assert audit['valid_transcripts']==2 and audit['invalid_transcripts']==0
    assert len(R.chainset(tmp_path/'B_predictions.gff3',{'plus':17,'minus':17}))==2


def test_single_cds_codon_coordinates():
    for strand in ('+','-'):
        parts=R.codon_pieces([(7,16)],strand)
        assert parts['start_codon']==([(7,10)] if strand=='+' else [(13,16)])
        assert parts['stop_codon']==([(13,16)] if strand=='+' else [(7,10)])
