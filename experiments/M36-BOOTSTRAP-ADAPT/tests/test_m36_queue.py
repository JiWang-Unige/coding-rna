"""M36 drivers on CPU: GPU gate (m36_gate.sh), the queue (m36_queue.sh) dry run on a FAKE root that stops before any GPU step, the CPU task queue engine (m36_cpu_queue.py),
the crash-resilient supervisor (m36_supervise.sh) and the retrying inference launcher (m36_infer.sh) with stub trainers / inferers.  No GPU, no model.
Run (ont): cd experiments/M36-BOOTSTRAP-ADAPT/tests && nice -n 10 taskset -c 25-31 ~/envs/generanno/bin/pytest -q -s test_m36_queue.py
"""
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import m36_common as W
import m36_cpu_queue as Q

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
ME = os.environ.get('USER') or subprocess.run(['id', '-un'], capture_output=True, text=True).stdout.strip()


def bash(cmd, env=None, timeout=300, cwd=None):
    return subprocess.run(['bash', '-c', cmd], capture_output=True, text=True, env=dict(os.environ, **(env or {})), timeout=timeout, cwd=cwd)


# ---- gate -------------------------------------------------------------------------------------------------------------------------------------------
@pytest.fixture
def procs():
    ps = []

    def start(argv0, cwd):
        p = subprocess.Popen(['bash', '-c', f'cd "{cwd}" && exec -a "{argv0}" sleep 120'])
        ps.append(p)
        time.sleep(0.4)
        return p
    yield start
    for p in ps:
        p.send_signal(signal.SIGTERM)
        p.wait()


def fake_smi(tmp_path, pids):
    d = tmp_path / 'bin'
    d.mkdir(exist_ok=True)
    (d / 'nvidia-smi').write_text('#!/bin/bash\n' + ''.join(f'echo {p}\n' for p in pids))
    (d / 'nvidia-smi').chmod(0o755)
    return str(d) + ':' + os.environ['PATH']


def blockers(tmp_path, pids=(), done=True, skip=False):
    f = tmp_path / 'P2B_DONE'
    if done:
        f.write_text('x')
    elif f.exists():
        f.unlink()
    env = dict(PATH=fake_smi(tmp_path, pids), P2B_DONE=str(f))
    if skip:
        env['M36_GATE_SKIP_PROCS'] = '1'
    r = bash(f'source {SCRIPTS}/m36_gate.sh; m36_blockers', env)
    assert r.returncode == 0, r.stderr
    return r.stdout.strip()


def test_gate_requires_p2b_done(tmp_path):
    assert 'no P2B_DONE' in blockers(tmp_path, done=False, skip=True)
    assert blockers(tmp_path, done=True, skip=True) == ''


def test_gate_blocks_on_other_coding_rna_jobs_only(tmp_path, procs):
    # (the real P2b / training processes of this machine may be running while the test runs: only the pids of the fake processes are asserted)
    p = procs('python m34_train_lora.py --arm T10', tmp_path)                       # a coding-rna training process
    r = blockers(tmp_path)
    assert 'processes:' in r and str(p.pid) in r
    assert blockers(tmp_path, skip=True) == ''                                      # dry-run knob: no process / GPU scan
    p.send_signal(signal.SIGTERM); p.wait()
    q = procs('python m34_p2b_infer.py --model x', tmp_path)                        # P2b itself still running: M36 must wait
    assert str(q.pid) in blockers(tmp_path)
    q.send_signal(signal.SIGTERM); q.wait()
    assert str(q.pid) not in blockers(tmp_path) and str(p.pid) not in blockers(tmp_path)


def test_gate_gpu_processes_only_count_when_they_belong_to_coding_rna(tmp_path, procs):
    mine = tmp_path / 'coding-rna-job'
    mine.mkdir()
    other = tmp_path / 'te-project'
    other.mkdir()
    a = procs('python job.py', mine)                                                # GPU process whose cwd is a coding-rna directory
    b = procs('python chrom_span_eval.py', other)                                   # the TE project's GPU process (same account): does not count
    r = blockers(tmp_path, pids=[a.pid, b.pid])
    assert f'gpu pids: {a.pid}' in r and str(b.pid) not in r
    assert 'gpu pids' not in blockers(tmp_path, pids=[b.pid])
    c = procs('python /home/x/coding-rna/outputs/run.py', other)                    # ... but a coding-rna command line does
    assert str(c.pid) in blockers(tmp_path, pids=[c.pid])
    r = bash(f'source {SCRIPTS}/m36_gate.sh; m36_gpu_blockers', dict(PATH=fake_smi(tmp_path, [a.pid, b.pid, 999999])))
    assert r.stdout.split() == [str(a.pid)]                                         # pids of other users / vanished pids never block


# ---- queue dry run -----------------------------------------------------------------------------------------------------------------------------------
def fake_roots(tmp_path, with_pot=True):
    o36, o34 = tmp_path / 'm36', tmp_path / 'm34'
    (o34 / 'pc2_T10').mkdir(parents=True)
    base = o34 / 'pc2_T10' / 'model_003072.pt'
    base.write_bytes(b'fake base model')
    (o34 / 'pc2_T10' / 'best_step.txt').write_text('3072\n')
    (o34 / 'p2b' / 'T10' / 'pot').mkdir(parents=True) if with_pot else (o34 / 'p2b' / 'T10').mkdir(parents=True)
    (o34 / 'p2b' / 'T10' / 'DONE_deadbeef.json').write_text('{}')
    h = hashlib.md5(base.read_bytes()).hexdigest()
    (o34 / 'p2b' / 'model_T10.json').write_text(json.dumps(dict(arm='T10', model=str(base), md5=h, bytes=15)))
    (o34 / 'p2b' / 'P2B_DONE').write_text('x')
    return o36, o34, h


def run_queue(tmp_path, o36, o34, extra_env=None, timeout=300):
    env = dict(M36_OUT=str(o36), M36_M34_OUT=str(o34), M36_LOGS=str(tmp_path / 'logs'), M36_DRY='1', M36_GATE_SKIP_PROCS='1', M36_SLEEP_OK='0.05', M36_SLEEP_WAIT='0.05',
               PATH=fake_smi(tmp_path, []))
    env.update(extra_env or {})
    return subprocess.run(['bash', str(SCRIPTS / 'm36_queue.sh')], env=dict(os.environ, **env), capture_output=True, text=True, timeout=timeout)


def events(o36):
    return [json.loads(x) for x in open(o36 / 'queue.jsonl')]


def test_queue_dry_run_walks_every_stage_without_running_anything(tmp_path):
    o36, o34, h = fake_roots(tmp_path)
    r = run_queue(tmp_path, o36, o34)
    assert r.returncode == 0, (r.stdout, r.stderr)
    ev = events(o36)
    names = [e['event'] for e in ev]
    assert names[:5] == ['queued', 'gate_open', 'base_verified', 'plan', 'would_run'] and names[-1] == 'queue_exit' and ev[0]['dry'] == 1
    plan = next(e for e in ev if e['event'] == 'plan')
    expected = ['S2:dryas_octopetala', 'S2:eschscholzia_californica', 'S1:dryas_octopetala', 'S3:dryas_octopetala', 'F1:dryas_octopetala', 'S1:eschscholzia_californica',
                'S3:eschscholzia_californica']
    assert plan['order'].split() == expected and plan['steps'] == 512
    would = [e for e in ev if e['event'] == 'would_run']
    stages = [(e['stage'], e.get('arm')) for e in would]
    assert stages[0] == ('cpu_seeds_queue', None) and stages[1][0] == 'zeroshot_A'
    assert [a for s, a in stages if s == 'train'] == [x.replace(':', '_') for x in expected]               # the seven trainings, in the planned order
    assert [a for s, a in stages if s == 'infer_B'] == [x.replace(':', '_') for x in expected]
    assert [s for s, _ in stages if s not in ('train', 'infer_B')] == ['cpu_seeds_queue', 'zeroshot_A', 'cpu_bpost_queue', 'score_main', 'score_post']
    z = next(e for e in would if e['stage'] == 'zeroshot_A')
    assert z['regions'] == 'NC_141366.1,NC_141368.1,NC_141370.1,NC_141372.1,NC_141374.1__val20M,NW_028238198.1,NW_028238420.1,NW_028238472.1__val20M'
    inf = {e['arm']: e['regions'] for e in would if e['stage'] == 'infer_B'}
    assert inf['S2_dryas_octopetala'] == 'NC_141367.1,NC_141369.1,NC_141371.1,NC_141373.1' and inf['S3_eschscholzia_californica'] == 'NW_028238309.1,NW_028238461.1,NW_028238483.1'
    assert 'train_start' not in names and 'zeroshot_start' not in names and 'infer_B_start' not in names and 'score_start' not in names
    assert json.load(open(o36 / 'base_model.json'))['md5'] == h
    assert not (o36 / 'zeroshot_A').exists() and not (o36 / 'train').exists() and not (o36 / 'evalB').exists() and not (o36 / 'M36_DONE').exists()
    # restart is idempotent: same base is accepted again
    assert run_queue(tmp_path, o36, o34).returncode == 0


def test_queue_regions_match_the_m30_split():
    text = (SCRIPTS / 'm36_queue.sh').read_text()
    import m30_common as M30
    A = [c for sp in W.SPECIES for c in M30.SPECIES[sp]['train'] + [M30.SPECIES[sp]['val_contig']]]
    assert f"A_REGIONS={','.join(A)}\n" in text
    for sp in W.SPECIES:
        assert f"[{sp}]={','.join(M30.SPECIES[sp]['B'])}" in text


def test_queue_refuses_a_different_base_model_and_a_missing_z_arm(tmp_path):
    o36, o34, h = fake_roots(tmp_path)
    (o34 / 'pc2_T10' / 'model_003072.pt').write_bytes(b'other model')
    r = run_queue(tmp_path, o36, o34)
    assert r.returncode == 11 and events(o36)[-1]['event'] == 'base_md5_mismatch'                           # not the model P2b used for T10 (= Z)
    o36, o34, h = fake_roots(tmp_path / 'b')
    (o34 / 'pc2_T10' / 'best_step.txt').write_text('2304\n')
    assert run_queue(tmp_path / 'b', o36, o34).returncode == 11 and events(o36)[-1]['event'] == 'base_step_mismatch'
    o36, o34, h = fake_roots(tmp_path / 'c')
    (o34 / 'p2b' / 'T10' / 'DONE_deadbeef.json').unlink()
    assert run_queue(tmp_path / 'c', o36, o34).returncode == 9 and events(o36)[-1]['event'] == 'z_arm_missing'
    o36, o34, h = fake_roots(tmp_path / 'd')
    o36.mkdir(parents=True, exist_ok=True)
    json.dump(dict(model='x', md5='0' * 32), open(o36 / 'base_model.json', 'w'))
    assert run_queue(tmp_path / 'd', o36, o34).returncode == 11 and events(o36)[-1]['event'] == 'base_changed'    # a restart with another base is refused


def test_queue_waits_for_p2b_done_and_is_single_instance(tmp_path):
    o36, o34, h = fake_roots(tmp_path)
    (o34 / 'p2b' / 'P2B_DONE').unlink()
    env = dict(os.environ, M36_OUT=str(o36), M36_M34_OUT=str(o34), M36_LOGS=str(tmp_path / 'logs'), M36_DRY='1', M36_GATE_SKIP_PROCS='1', M36_SLEEP_OK='0.05', M36_SLEEP_WAIT='0.2',
               PATH=fake_smi(tmp_path, []))
    p = subprocess.Popen(['bash', str(SCRIPTS / 'm36_queue.sh')], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        time.sleep(3)
        assert p.poll() is None                                                                           # still waiting: no P2B_DONE
        ev = events(o36)
        assert [e['event'] for e in ev] == ['queued', 'waiting'] and 'no P2B_DONE' in ev[1]['reasons']
        second = subprocess.run(['bash', str(SCRIPTS / 'm36_queue.sh')], env=env, capture_output=True, text=True, timeout=60)
        assert second.returncode == 3 and 'already running' in second.stdout
        (o34 / 'p2b' / 'P2B_DONE').write_text('x')
        out, err = p.communicate(timeout=120)
        assert p.returncode == 0, (out, err)
        names = [e['event'] for e in events(o36)]
        assert names[:3] == ['queued', 'waiting', 'gate_open'] and names[-1] == 'queue_exit'
    finally:
        if p.poll() is None:
            p.kill()


# ---- CPU task queue ------------------------------------------------------------------------------------------------------------------------------------
def make_task(tmp, name, seq, strand, length, group, fail=False, delay=0.):
    run = tmp / 'run'
    out = tmp / 'out' / name
    done = out / f'{seq}_{strand}.json'
    cmd = ['bash', '-c', f'sleep {delay}; ' + ('exit 3' if fail else f'mkdir -p {out}; echo {{}} > {done}')]
    return Q.Task('post', group, run, seq, strand, out, length, cmd, done)


def test_cpu_queue_waits_for_inputs_largest_first_and_fires_group_hooks(tmp_path):
    ev = []
    run = tmp_path / 'run'
    (run / 'pot').mkdir(parents=True)
    tasks = [make_task(tmp_path, 'g1', 'cA', '+', 10, 'g1'), make_task(tmp_path, 'g1', 'cB', '+', 30, 'g1'), make_task(tmp_path, 'g2', 'cC', '+', 20, 'g2')]
    fired = []
    hooks = {'g1': lambda: fired.append('g1'), 'g2': (lambda: fired.append('g2'), lambda: True)}                # g2's hook is skipped ("already done")

    def deliver():
        time.sleep(0.6)
        for seq, L in (('cA', 10), ('cB', 30), ('cC', 20)):                                                    # pot first, chains later (the decoder writes chains last)
            (run / 'pot' / f'{seq}_+.npy').write_bytes(b'x')
        time.sleep(0.4)
        (run / 'chains_tmp').mkdir()
        for seq in ('cA', 'cB', 'cC'):
            (run / 'chains_tmp' / f'{seq}_+.json').write_text('{"chains": []}')
        os.rename(run / 'chains_tmp', run / 'chains')                                                              # all three decoded chains appear at once
    th = threading.Thread(target=deliver)
    th.start()
    failed, done_hooks = Q.run_tasks(tasks, hooks, workers=1, mem_budget=10**12, cpus=None, emit=lambda **k: ev.append(k), poll=0.05)
    th.join()
    assert failed == [] and fired == ['g1'] and done_hooks == ['g1', 'g2']
    launches = [e['task'] for e in ev if e['event'] == 'launch']
    assert launches == ['post:g1:cB:+', 'post:g2:cC:+', 'post:g1:cA:+']                                          # nothing starts before its pot AND chains exist; then largest first
    assert all(not t.ready() or t.done() for t in tasks)


def test_cpu_queue_memory_budget_and_worker_cap(tmp_path):
    run = tmp_path / 'run'
    (run / 'pot').mkdir(parents=True)
    (run / 'chains').mkdir()
    tasks = [make_task(tmp_path, 'g', f'c{i}', '+', 100, 'g', delay=0.4) for i in range(4)]
    for t in tasks:
        (run / 'pot' / f'{t.seq}_+.npy').write_bytes(b'x')
        (run / 'chains' / f'{t.seq}_+.json').write_text('{}')
    ev = []
    t0 = time.time()
    Q.run_tasks(tasks, {}, workers=4, mem_budget=Q.MEM_PER_BASE['post'] * 250, cpus=None, emit=lambda **k: ev.append((time.time() - t0, k)), poll=0.05)   # budget for two 100-base tasks
    starts = [t for t, k in ev if k['event'] == 'launch']
    assert len(starts) == 4 and starts[2] - starts[0] >= 0.35 and starts[1] - starts[0] < 0.3                      # two at a time


def test_cpu_queue_retries_once_then_reports_the_failure(tmp_path):
    run = tmp_path / 'run'
    (run / 'pot').mkdir(parents=True)
    (run / 'chains').mkdir()
    bad = make_task(tmp_path, 'g', 'bad', '+', 5, 'g', fail=True)
    good = make_task(tmp_path, 'g', 'ok', '+', 4, 'g')
    for t in (bad, good):
        (run / 'pot' / f'{t.seq}_+.npy').write_bytes(b'x')
        (run / 'chains' / f'{t.seq}_+.json').write_text('{}')
    ev = []
    failed, _ = Q.run_tasks([bad, good], {'g': lambda: ev.append(('hook',))}, workers=2, mem_budget=10**12, cpus=None, emit=lambda **k: ev.append(k), poll=0.05)
    assert failed == ['post:g:bad:+'] and good.done() and [e['event'] for e in ev if isinstance(e, dict) and e['event'] == 'fail'] == ['fail', 'fail']
    assert ('hook',) not in ev                                                                                     # the group is not complete: its hook (seed build) does not fire


def test_cpu_queue_task_commands_use_the_real_scripts(tmp_path):
    t = Q.post_task('dryas_octopetala', tmp_path / 'run', 'NC_141366.1', '+', 1.0, tmp_path / 'post', 33_681_925)
    assert t.cmd[1:] == [str(Q.M33_POST), str(tmp_path / 'run'), 'NC_141366.1', '+', '1', str(tmp_path / 'post')] and t.done_file == tmp_path / 'post' / 'tau1' / 'NC_141366.1_+.json'
    n = Q.nogene_task('dryas_octopetala', tmp_path / 'run', 'NC_141366.1', '-', 1.0, tmp_path / 'ng', 33_681_925)
    assert n.cmd[1:] == [str(SCRIPTS / 'm36_nogene_marginal.py'), 'run', str(tmp_path / 'run'), 'NC_141366.1', '-', '1', str(tmp_path / 'ng')]
    assert n.done_file == tmp_path / 'ng' / 'tau1' / 'nogene_NC_141366.1_-.json' and Q.M33_POST.exists() and (SCRIPTS / 'm36_nogene_marginal.py').exists()
    assert t.mem() < n.mem()
    tasks = Q.seeds_tasks(tmp_path / 'r', tmp_path / 'p', tmp_path / 'n', list(W.SPECIES))
    assert len(tasks) == 2 * 2 * (5 + 3) and {x.seq for x in tasks} == {c for sp in W.SPECIES for c in W.a_lengths(sp)}      # (post + nogene) x 2 strands x 8 A contigs
    import m30_common as M30
    assert not {x.seq for x in tasks} & {c for sp in W.SPECIES for c in M30.SPECIES[sp]['B']}                       # B contigs never appear
    bt, hooks = Q.bpost_tasks([('S2_dryas_octopetala', str(tmp_path / 'e'), 'dryas_octopetala', False), ('Z_dryas_octopetala', str(tmp_path / 'z'), 'dryas_octopetala', True)], tmp_path / 'postB')
    assert len(bt) == 2 * 4 * 2 and set(hooks) == {'S2_dryas_octopetala'} and {x.seq for x in bt} == set(M30.SPECIES['dryas_octopetala']['B'])


def test_bpost_hook_deletes_pot_but_never_the_z_potentials(tmp_path):
    for name, keep in (('arm', False), ('z', True)):
        (tmp_path / name / 'pot').mkdir(parents=True)
        (tmp_path / name / 'pot' / 'x.npy').write_bytes(b'x')
    tasks, hooks = Q.bpost_tasks([('arm_dryas_octopetala', tmp_path / 'arm', 'dryas_octopetala', False), ('z_dryas_octopetala', tmp_path / 'z', 'dryas_octopetala', True)], tmp_path / 'postB')
    hooks['arm_dryas_octopetala']()
    assert not (tmp_path / 'arm' / 'pot').exists() and (tmp_path / 'z' / 'pot' / 'x.npy').exists()


# ---- supervisor / inference launcher with stubs --------------------------------------------------------------------------------------------------------
STUB_TRAIN = r'''
import json, sys, os
import torch
a = sys.argv[1:]
out = a[a.index('--out') + 1]; stop = int(a[a.index('--stop-at') + 1])
state = os.path.join(out, 'stub_attempts')
n = int(open(state).read()) + 1 if os.path.exists(state) else 1
open(state, 'w').write(str(n))
mode = os.environ['STUB_MODE']
if mode == 'crash_once' and n == 1:
    torch.save(dict(step=128), os.path.join(out, 'ckpt_latest.pt')); sys.exit(1)
if mode == 'always_fail':
    sys.exit(2)
if mode == 'refuse':
    sys.exit(11)
ck = os.path.join(out, 'ckpt_latest.pt')
start = torch.load(ck)['step'] if os.path.exists(ck) else 0
json.dump(dict(completed=True, start=start, args=a), open(os.path.join(out, f'DONE_{stop}.json'), 'w'))
'''


def run_supervise(tmp_path, mode, stop=512):
    tmp_path.mkdir(parents=True, exist_ok=True)
    stub = tmp_path / 'stub_train.py'
    stub.write_text(STUB_TRAIN)
    out = tmp_path / 'train' / 'S2_x'
    env = dict(STUB_MODE=mode, M36_PY_TRAIN=sys.executable, M36_TRAIN_SCRIPT=str(stub), M36_FAKE_FREE_MIB='20000', M36_LOGS=str(tmp_path / 'logs'), M36_SCRIPTS=str(tmp_path),
               M34_CPUS='25-31', M36_SLEEP_RETRY='0', M36_SLEEP_BUSY='0')
    r = bash(f'bash {SCRIPTS}/m36_supervise.sh S2 x {out} {stop}', env)
    ev = [json.loads(x) for x in open(out / 'interruptions.jsonl')]
    return r, ev, out


def test_supervisor_resumes_after_a_crash(tmp_path):
    r, ev, out = run_supervise(tmp_path, 'crash_once')
    assert r.returncode == 0, r.stderr
    assert [e['event'] for e in ev] == ['launch', 'interruption', 'launch', 'done']
    assert ev[1]['last_checkpoint_step'] == '128' and ev[1]['rc'] == 1
    d = json.load(open(out / 'DONE_512.json'))
    assert d['start'] == 128 and '--resume' in d['args'] and d['args'][d['args'].index('--stop-at') + 1] == '512'                  # the second attempt resumed from the checkpoint
    assert ['--arm', 'S2', '--species', 'x'] == d['args'][d['args'].index('--arm'):d['args'].index('--arm') + 4]


def test_supervisor_gives_up_on_a_deterministic_error_and_passes_a_refusal_through(tmp_path):
    r, ev, out = run_supervise(tmp_path, 'always_fail')
    assert r.returncode == 9 and ev[-1]['event'] == 'FATAL' and 3 <= sum(e['event'] == 'launch' for e in ev) <= 4   # three failed attempts without progress -> FATAL (the M34 supervisor's rule)
    r, ev, out = run_supervise(tmp_path / 'b', 'refuse')
    assert r.returncode == 11 and [e['event'] for e in ev] == ['launch', 'refused']


def test_supervisor_waits_while_the_gpu_is_busy(tmp_path):
    stub = tmp_path / 'stub_train.py'
    stub.write_text(STUB_TRAIN)
    out = tmp_path / 'train' / 'S2_x'
    env = dict(STUB_MODE='crash_once', M36_PY_TRAIN=sys.executable, M36_TRAIN_SCRIPT=str(stub), M36_FAKE_FREE_MIB='5000', M36_LOGS=str(tmp_path / 'logs'), M36_SCRIPTS=str(tmp_path),
               M34_CPUS='25-31', M36_SLEEP_RETRY='0', M36_SLEEP_BUSY='0')
    r = bash(f'bash {SCRIPTS}/m36_supervise.sh S2 x {out} 512', env, timeout=120)                                         # < 7.5 GB free on every attempt: never launches
    ev = [json.loads(x) for x in open(out / 'interruptions.jsonl')]
    assert r.returncode == 9 and all(e['event'] in ('gpu_busy_wait', 'FATAL') for e in ev) and not (out / 'stub_attempts').exists()


STUB_INFER = r'''
import json, sys, os
a = sys.argv[1:]
out = a[a.index('--out') + 1]
os.makedirs(out, exist_ok=True)
state = os.path.join(out, 'stub_attempts')
n = int(open(state).read()) + 1 if os.path.exists(state) else 1
open(state, 'w').write(str(n))
if n == 1:
    sys.exit(1)
json.dump(dict(completed=True, argv=a), open(os.path.join(out, 'DONE_abcd.json'), 'w'))
'''


def test_infer_launcher_retries_and_keeps_the_potentials(tmp_path):
    stub = tmp_path / 'stub_infer.py'
    stub.write_text(STUB_INFER)
    out = tmp_path / 'evalB' / 'S2_x'
    env = dict(M36_PY_TRAIN=sys.executable, M36_INFER_SCRIPT=str(stub), M36_FAKE_FREE_MIB='20000', M36_LOGS=str(tmp_path / 'logs'), M36_M34_SCRIPTS=str(tmp_path), M34_CPUS='25-31',
               M36_SLEEP_RETRY='0', M36_SLEEP_BUSY='0')
    r = bash(f'bash {SCRIPTS}/m36_infer.sh {out} /m/model.pt NC_1,NC_2 4', env)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out / 'DONE_abcd.json'))['argv']
    assert d[d.index('--model') + 1] == '/m/model.pt' and d[d.index('--regions') + 1] == 'NC_1,NC_2' and '--keep-head' in d and d[d.index('--decode-workers') + 1] == '4'
    ev = [json.loads(x) for x in open(out / 'interruptions.jsonl')]
    assert [e['event'] for e in ev] == ['interruption', 'done']
    r = bash(f'bash {SCRIPTS}/m36_infer.sh {tmp_path}/o2 /m/model.pt NC_1 4', dict(env, M36_KEEP_POT='0'))
    assert '--keep-head' not in json.load(open(tmp_path / 'o2' / 'DONE_abcd.json'))['argv']
