import importlib.util
from pathlib import Path
import numpy as np
import pytest

BASE=Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY"


def module(name):
    spec=importlib.util.spec_from_file_location(name,BASE/(name+".py"))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T=module("wt_tiberius")
C=module("wt_compare")


def test_cache_replay_keeps_scores_and_never_calls_network(tmp_path):
    x=np.eye(5,dtype=np.float32)[np.array([[0,1,2]])]
    y=np.arange(45,dtype=np.float32).reshape(1,3,15)
    native=T.NeuralCache(tmp_path/"cache","record")
    native.call(x,{"strand":"+"},lambda:y)
    native.finish()
    replay=T.NeuralCache(tmp_path/"cache","replay")
    def forbidden():
        raise AssertionError("Replay must not infer")
    assert np.array_equal(replay.call(x,{"strand":"+"},forbidden),y)
    replay.finish()


def test_changed_input_and_missing_request_fail(tmp_path):
    x=np.zeros((1,3,5),dtype=np.float32)
    cache=T.NeuralCache(tmp_path/"cache","record")
    cache.call(x,{"strand":"+"},lambda:np.zeros((1,3,15),dtype=np.float32))
    cache.finish()
    replay=T.NeuralCache(tmp_path/"cache","replay")
    with pytest.raises(RuntimeError,match="input differs"):
        replay.call(x+1,{"strand":"+"},lambda:None)
    replay=T.NeuralCache(tmp_path/"cache","replay")
    with pytest.raises(RuntimeError,match="fewer requests"):
        replay.finish()
    replay.call(x,{"strand":"+"},lambda:None)
    with pytest.raises(RuntimeError,match="uncached"):
        replay.call(x,{"strand":"+"},lambda:None)


def test_cache_budget_stops_before_writing_scores(tmp_path):
    cache=T.NeuralCache(tmp_path/"cache","record",limit=1)
    with pytest.raises(RuntimeError,match="limit"):
        cache.call(np.zeros((1,2,5)),{},lambda:np.zeros((1,2,15)))
    assert not list((tmp_path/"cache").glob("*.npy"))


def test_exact_matching_phase_strand_and_duplicate_policy():
    row={"strand":"-","CDS_chain_0based_halfopen":[[10,15],[30,34]]}
    assert C.reference_key(row)==("chr22","-",((10,15,2),(30,34,0)))
    tx={"id":"x","seqid":"chr22","strand":"-","CDS":[(30,34,"0"),(10,15,"2")]}
    assert C.match_target(row,[tx])["unique_WT_exact"]
    assert not C.match_target(row,[{**tx,"strand":"+"}])["unique_WT_exact"]
    assert C.match_target(row,[tx,{**tx,"id":"y"}])["ambiguous_duplicate_exact"]


def test_prediction_parser_preserves_native_CDS_with_terminal_stop(tmp_path):
    path=tmp_path/"a.gtf"
    path.write_text('chr22\tTiberius\tCDS\t11\t15\t.\t-\t2\tgene_id "g"; transcript_id "t";\n'
                    'chr22\tTiberius\tCDS\t31\t34\t.\t-\t0\tgene_id "g"; transcript_id "t";\n'
                    'chr22\tTiberius\tstop_codon\t11\t13\t.\t-\t0\tgene_id "g"; transcript_id "t";\n')
    tx=C.parsed(path)
    assert len(tx)==1 and C.key(tx[0])==("chr22","-",((10,15,2),(30,34,0)))
