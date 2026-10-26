"""CPU-only smoke of the REAL `m36_train.main` -> `m34_train_lora.main` loop (no GPU, no real encoder), used by tests/test_m36_train.py.  Run under the pc2 backend:
  M34_BACKEND=pc2 ~/envs/plantcad2/bin/python m36_cpu_smoke.py <workdir> <arm> <stop_at> [resume]
Everything that would touch CUDA or the PlantCAD2 encoder is replaced by a CPU fake:
  * `torch.Tensor.to / Module.to('cuda')` -> no-op; `torch.cuda.synchronize / max_memory_*` -> 0;
  * `BK.load` -> a fake ctx whose 'encoder' is a small module of trainable 'LoRA' parameters (+ a fake `lora_api`); `m36_train.Patched` wraps it and must load the BASE model's LoRA into it;
  * `BK._pc2_encode` -> RC-equivariant fake features (as tests/test_m34_pc2_backend.py), made to depend on the LoRA parameters so that gradients reach them;
  * `BK.window_pots` -> the real function with device='cpu'.
The target species is a synthetic 2.4-kb contig with planted genes; `m36_train.main` gets it through its `hooks` (inputs / labeler); the draw schedule, the seed labeler (all four arms),
`T.example`, the head, the parallel window DP (real process pool), AdamW + cosine LR, checkpoint save / resume and the model files are the production code.
Writes <workdir>/{data, train_<arm>/...} and <workdir>/first_lora_<arm>_<phase>.json (LoRA parameters seen by the first forward pass of the run).
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T      # noqa: E402  (first: backend adapter)
import numpy as np                # noqa: E402
import torch                      # noqa: E402

assert T.BACKEND == 'pc2', 'run with M34_BACKEND=pc2 under ~/envs/plantcad2'
import m34_backend as BK          # noqa: E402
import m34_train_lora as TL       # noqa: E402
import m36_common as W            # noqa: E402
import m36_seed_labels as SL      # noqa: E402
import m36_schedule as SCH        # noqa: E402
import m36_train as TR            # noqa: E402
import manifest as MF             # noqa: E402

SPECIES = 'dryas_octopetala'      # (the synthetic contig is short: the head is a BiGRU over the bases, ~10 s per 1 kb on a loaded CPU)
SAFE = ['GCT', 'GAA', 'AAA', 'CCT', 'GGC', 'TTC', 'GAC', 'CAA']
FAKE_LORA_CFG = dict(rank=2, fake=True)


def make_gene(rng, n_exons):
    n_codons = int(rng.integers(12, 30))
    cds = 'ATG' + ''.join(rng.choice(SAFE, n_codons)) + 'TAA'
    cuts = sorted(rng.choice(np.arange(4, len(cds) - 4), n_exons - 1, replace=False)) if n_exons > 1 else []
    pieces, prev = [], 0
    for c in list(cuts) + [len(cds)]:
        pieces.append(cds[prev:c]); prev = c
    dna, exons, pos = '', [], 0
    for j, piece in enumerate(pieces):
        exons.append((pos, pos + len(piece))); dna += piece; pos += len(piece)
        if j < len(pieces) - 1:
            intron = 'GT' + ''.join(rng.choice(list('ACGT'), int(rng.integers(8, 30)))) + 'AG'
            dna += intron; pos += len(intron)
    return dna, exons


def synthetic(workdir, length=2400, n_genes=3, seed=7):
    rng = np.random.default_rng(seed)
    dna = list(rng.choice(list('ACGT'), length))
    rows, cursor = [], 200
    for k in range(n_genes):
        g, exons = make_gene(rng, int(rng.integers(1, 4)))
        strand = '+' if k % 2 == 0 else '-'
        if strand == '+':
            dna[cursor:cursor + len(g)] = list(g)
            cds = [(cursor + a, cursor + b) for a, b in exons]
        else:
            dna[cursor:cursor + len(g)] = list(MF.rc(g))
            cds = sorted((cursor + len(g) - b, cursor + len(g) - a) for a, b in exons)
        rows.append(SL.seed_row('SYN', strand, SL.phases(cds, strand), f'syn{k}'))
        cursor += len(g) + 250
    return ''.join(dna), rows


_TABLE = {}


def fake_features(strings, device='cpu'):
    """RC-equivariant fake per-base features [len, 1536] bf16 (v(comp b) = reverse(v(b)); N -> zeros), built by table lookup (fast enough for the CPU smoke)"""
    if not _TABLE:
        rng = np.random.default_rng(7)
        X = {b: torch.from_numpy(rng.standard_normal(768).astype(np.float32)) for b in 'ACGT'}
        comp = dict(zip('ACGT', 'TGCA'))
        rows = [torch.cat([X[c], X[comp[c]].flip(0)]) for c in 'ACGT'] + [torch.zeros(1536)]
        _TABLE['t'] = torch.stack(rows).to(torch.bfloat16)
        lut = np.full(256, 4, np.int64)
        for i, c in enumerate('ACGT'):
            lut[ord(c)] = i
        _TABLE['lut'] = lut
    return [_TABLE['t'][torch.from_numpy(_TABLE['lut'][np.frombuffer(s.encode(), dtype=np.uint8)])] for s in strings]


class FakeEnc(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.lora_A = torch.nn.Parameter(torch.zeros(4))
        self.lora_B = torch.nn.Parameter(torch.zeros(3))


class FakeLoraApi:
    CFG = FAKE_LORA_CFG

    @staticmethod
    def lora_state(enc):
        return {k: p.detach().cpu().clone() for k, p in enc.named_parameters()}

    @staticmethod
    def load_lora_state(enc, state, cfg):
        with torch.no_grad():
            for k, p in enc.named_parameters():
                p.copy_(state[k])


def install_fakes(workdir, tag):
    _to, _mto = torch.Tensor.to, torch.nn.Module.to
    cpu = lambda a, k: (tuple('cpu' if x == 'cuda' else x for x in a), {kk: ('cpu' if vv == 'cuda' else vv) for kk, vv in k.items()})
    torch.Tensor.to = lambda self, *a, **k: _to(self, *cpu(a, k)[0], **cpu(a, k)[1])
    torch.nn.Module.to = lambda self, *a, **k: _mto(self, *cpu(a, k)[0], **cpu(a, k)[1])
    torch.cuda.synchronize = lambda *a, **k: None
    torch.cuda.max_memory_allocated = lambda *a, **k: 0
    torch.cuda.max_memory_reserved = lambda *a, **k: 0
    real_window_pots = BK.window_pots
    seen = {}

    def fake_load(cfg, lora=True, device='cuda'):
        enc = FakeEnc()
        cfg['lora'] = dict(FAKE_LORA_CFG)
        return BK.Ctx(enc, cfg, FakeLoraApi), list(enc.parameters())

    def fake_window_pots(head, ctx, it, idx, grad, device='cuda'):
        if 'first' not in seen:                                                # LoRA values seen by the first forward pass of this run
            seen['first'] = {k: v.tolist() for k, v in FakeLoraApi.lora_state(ctx.enc).items()}
            json.dump(seen['first'], open(Path(workdir) / f'first_lora_{tag}.json', 'w'))
        pots = real_window_pots(head, ctx, it, idx, grad, device='cpu')
        scale = 1. + 1e-2 * sum(p.sum() for p in ctx.enc.parameters())          # the 'encoder' depends on its LoRA parameters -> gradients reach them
        return [p * (scale if grad else scale.detach()) for p in pots]

    BK.load = fake_load
    BK.window_pots = fake_window_pots
    BK._pc2_encode = lambda enc, strings, layer=-1, batch=8, grad=False, device='cpu': fake_features(strings)


def prepare(workdir):
    """base model file (head + non-zero fake LoRA), synthetic genome / seeds -> dict of paths"""
    d = Path(workdir) / 'data'
    d.mkdir(parents=True, exist_ok=True)
    base = d / 'base_model.pt'
    dna, rows = synthetic(workdir)
    fa = d / 'genome.fa'
    if not fa.exists():
        fa.write_text('>SYN\n' + '\n'.join(dna[i:i + 80] for i in range(0, len(dna), 80)) + '\n')
    if not base.exists():
        torch.manual_seed(0)
        head = BK.new_head(dict(layers=None), 0)
        cfg = dict(BK.encoder_cfg(dict(layer=-1)), lora=dict(FAKE_LORA_CFG))
        lora = dict(lora_A=torch.tensor([.5, -.25, .125, 1.]), lora_B=torch.tensor([.3, -.6, .9]))
        torch.save(dict(model=head.state_dict(), lora=lora, step=3072, config=cfg), base)
    return dict(base=base, fa=fa, rows=rows, lens={'SYN': len(dna)})


def main():
    workdir, arm, stop, resume = sys.argv[1], sys.argv[2], int(sys.argv[3]), len(sys.argv) > 4 and sys.argv[4] == 'resume'
    P = prepare(workdir)
    os.makedirs(workdir, exist_ok=True)
    sched = Path(workdir) / 'data' / f'schedule_{arm}.jsonl'
    if not sched.exists():
        SCH.write_schedule(SCH.draw_schedule(SPECIES, P['rows'], P['lens'], steps=3, batch=1), sched)
    install_fakes(workdir, f"{arm}_{'resume' if resume else 'start'}")
    rows = P['rows']
    bg = {('SYN', '+'): [(0, 150)], ('SYN', '-'): [(P['lens']['SYN'] - 150, P['lens']['SYN'])]} if arm in ('S1', 'S3') else None
    hooks = dict(inputs=lambda a, sp: dict(arm=a, species=sp, rows=rows, s1=None, s2=None, background=bg, dropped_s1_overlapping_s2=0, files={}),
                 labeler=lambda inputs: SL.SeedLabeler(SPECIES, P['fa'], P['lens'], inputs['rows'], inputs['background']))
    argv = ['run', '--arm', arm, '--species', SPECIES, '--out', str(Path(workdir) / f'train_{arm}'), '--base', str(P['base']), '--schedule', str(sched),
            '--steps', '3', '--batch', '1', '--stop-at', str(stop), '--workers', '2', '--ckpt-every', '1', '--keep-every', '2'] + (['--resume'] if resume else [])
    TR.main(argv, hooks)


if __name__ == '__main__':
    main()
