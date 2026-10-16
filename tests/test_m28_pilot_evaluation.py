import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_pilot as I
import evaluate_pilot as V

def test_inference_grid_rejects_missing_duplicate_or_out_of_scope_windows():
    data=[]
    for species,lengths in I.DEV_LENGTHS.items():
        for q,length in lengths.items():
            for strand in ('+','-'):
                for a in I.M.starts(length):
                    data.append({'split':'val','window_id':species+'|'+q+'|'+str(a)+'|'+strand})
    groups=I.dev_groups(data)
    assert len(groups)==4 and sum(map(len,groups.values()))==8690
    for invalid in (data[:-1],data[:-1]+data[:1],data+[{'split':'val','window_id':'sealed|q|0|+'}]):
        with pytest.raises(ValueError): I.dev_groups(invalid)

def test_reference_attrition_uses_actual_chain_keys_not_IDs_or_injected_positives():
    refs=[{'CDS':[(0,9,'0')]},{'CDS':[(20,29,'0')]},{'CDS':[(40,49,'0')]}]
    assembled={'stages':{'free_unique':[[[0,9]],[[20,29]]],'owner_unique':[[[0,9]]],'selected_unique':[]},
               'counts':{'free_records':2}}
    audit=V.stage_audit(assembled,refs)
    assert audit['not_in_actual_free_candidates']==1
    assert audit['lost_at_owner']==1
    assert audit['lost_at_conflict_or_nonpositive_gain']==1

def test_macro_is_per_species_mean_not_micro_or_harmonic_of_macro_PR():
    def result(tp,pred,ref):
        metric=V.E.S.prf(tp,pred,ref)
        return {'exact_chain':metric,'CDS_base_unstranded_union':metric,'historical_parent_exact_chain':metric,
                'background':{'bp':10,'predicted_span_bp':1,'predicted_CDS_bp':1,'wholly_background_unique_chains':1}}
    a,b=result(1,1,1),result(0,9,9)
    pooled=V.aggregate({'a':a,'b':b})
    assert pooled['exact_chain']['f1']==.1
    assert pooled['macro_exact_chain']['f1']==.5
