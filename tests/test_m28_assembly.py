import sys
import itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.assembly import assemble_predictions,nearest_owner
from src.m28.c0_decode import owner_window
from src.m28.core import select_nonoverlapping

def test_tail_ownership_matches_existing_rule():
    windows=[(0,24),(12,36),(17,41)]
    for a in range(41):
        for b in range(a+1,42):
            assert nearest_owner([(a,b)],windows)==owner_window([(a,b)],windows)

def test_global_stage_ledger_nonowner_loss_duplicates_and_conflict():
    w=[(0,24),(12,36)]
    records=[{'window_index':0,'chains':[[(2,8)],[(2,8)],[(16,24)]],'score':[2,1,5]},
             {'window_index':1,'chains':[[(8,18)],[(10,20)]],'score':[3,4]}]
    r=assemble_predictions(records,w,'+')
    assert r['counts']=={'free_records':5,'free_unique':4,'owner_records':4,'owner_unique':3,
                        'duplicate_owner_records':1,'nonowner_only_unique_lost':1,'selected_unique':2,'score_field':'score'}
    assert [x['chain'] for x in r['predictions']]==[((2,8),),((22,32),)]
    minus=assemble_predictions([{'window_index':0,'chains':[[(16,22)]],'score':[2]}],w,'-')
    assert minus['predictions'][0]['chain']==((2,8),)

def test_linear_memory_selector_matches_exhaustive_positive_gain():
    chains=[[(0,6)],[(4,12)],[(6,9)],[(12,18)],[(20,24)],[(8,21)]]
    for scores in itertools.product((-1,0,2),repeat=len(chains)):
        selected=select_nonoverlapping(chains,scores)
        best=0
        for mask in itertools.product((0,1),repeat=len(chains)):
            ids=[i for i,v in enumerate(mask) if v]
            if all(chains[i][-1][1]<=chains[j][0][0] or chains[j][-1][1]<=chains[i][0][0] for i,j in itertools.combinations(ids,2)):
                best=max(best,sum(scores[i] for i in ids))
        assert sum(scores[i] for i in selected)==best
