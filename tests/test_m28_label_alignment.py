import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from align_labels import relabel_window

def test_updated_primary_changes_target_not_geometry_or_transcript_protection():
    w={'start':0,'end':1000,'strand':'-','positive_ids':['old'],'overlapping_transcript_ids':['old','new'],
       'chain_unknown_intervals':[],'fine_label_unknown_intervals':[]}
    c={'old':{'id':'old','span':[300,600],'eligible_positive':False,'CDS_complete_supported':False,'overlapping_primary':False},
       'new':{'id':'new','span':[350,620],'eligible_positive':True,'CDS_complete_supported':True,'overlapping_primary':False}}
    changed=relabel_window(w,c)
    assert changed['positive_ids']==['new']
    assert changed['fine_label_unknown_intervals']==[[400,700]]
    assert changed['chain_unknown_intervals']==[[400,700]]
    assert changed['overlapping_transcript_ids']==w['overlapping_transcript_ids']
    assert w['positive_ids']==['old']
