import json
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION'))
import run_graph_run as R


@pytest.mark.parametrize('strand',['+','-'])
@pytest.mark.parametrize('empty',[False,True])
def test_json_roundtrip_uses_frozen_metric(tmp_path,strand,empty):
    key=('synthetic','chr1')
    model={'strand':strand,'cds':[(3,6),(12,18)],'phase':[0,0]}
    predictions={key:[] if empty else [model]}
    references={key:[{'strand':strand,'CDS':[(3,6,'0'),(12,18,'0')],'exon':[]}]}
    lengths={key:30}
    expected=R.D.m25._validation_metrics(predictions,references,lengths)
    path=tmp_path/'models.json'
    path.write_text(json.dumps([{'species':key[0],'seqid':key[1],'models':predictions[key]}]))
    loaded=R.load_predictions(path)
    assert loaded==predictions
    assert R.D.m25._validation_metrics(loaded,references,lengths)==expected
    assert expected['exact_CDS_chain_F1']==(0.0 if empty else 1.0)
    assert expected['gene_count_ratio']==(0.0 if empty else 1.0)


def test_saved_score_source_is_read_only(tmp_path,monkeypatch):
    source=tmp_path/'original'
    (source/'scores').mkdir(parents=True)
    monkeypatch.setattr(R,'OUT',tmp_path/'recovery')
    monkeypatch.setattr(R,'SCORE_SOURCE',source)
    monkeypatch.setattr(R,'SCOPE',(('synthetic','chr1',12),))
    for name,width in (('region',3),('boundary',4),('phase',4)):
        np.save(R.score_path('chr1','+',name,source),np.ones((12,width),dtype=np.float16))
    arrays=R.load_scores('chr1','+')
    assert [a.shape for a in arrays]==[(12,3),(12,4),(12,4)]
    assert all(not a.flags.writeable for a in arrays)
    assert not R.OUT.exists()
