"""Switchable encoder backend of the M34 trainer / inference (env M34_BACKEND = ntv2 | pc2; chosen before `m34_train_common` is imported).

  ntv2 : NT-v2-500M-multi-species (6-mer tokens, 12,288-bp windows, stride 6,144), hidden_states 21+27 concatenated, q/v LoRA on blocks 13..27 -- the M31 B recipe
         (code: M31 `ntv2_adapter` / `c1_model`, unchanged).  Env: ~/envs/generanno.
  pc2  : PlantCAD2-Small (single-nucleotide, 8,192-bp windows, stride 6,144, ONE forward gives both strands by exact RC equivariance, last_hidden_state,
         autocast bf16, per-base features mean-pooled to six-base tokens before the frozen C0 head of width 8,196) -- the M32 geometry (report 17 section 1).
         Code: M32 `pc2_adapter` / `c1_model`, unchanged.  Env: ~/envs/plantcad2.
         LoRA on PlantCAD2 is NOT defined by this file: it is taken from the M35 experiment (contract below) so that M34 uses exactly the adaptation M35 validated.
         Contract: experiments/M35-PLANTCAD2-LORA/scripts/pc2_lora.py exposes  add_lora(model, cfg) -> [trainable params],  lora_state(model) -> {key: tensor},
         load_lora_state(model, state, cfg),  CFG (dict, JSON-able; stored in every checkpoint).  Without it (`--lora none`) the PlantCAD2 encoder stays frozen and only the head trains
         (= M32 recipe on the M34 multi-species schedule).

The trainer / inference call only the functions below.  Draw items carry `seqs` (oriented window strings, the head's DNA input) and `gen` (genomic strings, PlantCAD2 only).
"""
import re

import torch

import m34_train_common as T

A, MD, G, Q, M = T.A, T.MD, T.G, T.Q, T.M
NT = T.BACKEND == 'ntv2'


# ---- configuration / construction -----------------------------------------------------------------------------------------------------------
def encoder_cfg(init_cfg):
    """Entries describing the encoder + head interface that go into the training config (and are checked at resume)."""
    if NT:
        return dict(backend='ntv2', encoder=A.MODEL, nt_width=A.NT_WIDTH, nt_stride=A.NT_STRIDE, layers=list(init_cfg['layers']),
                    layer_norm_per_layer=init_cfg['layer_norm_per_layer'], lora_r=8, lora_alpha=16, lora_encoder_layers=list(range(12, 27)))
    return dict(backend='pc2', encoder=A.MODEL, pc_width=A.PC_WIDTH, pc_stride=A.PC_STRIDE, head_width=A.HEAD_W, layer=init_cfg['layer'],
                token_pool='mean over valid bases of aligned 6-base groups')


def check_init(init):
    if NT:
        assert init['step'] == 2048 and tuple(init['config']['layers']) == (21, 27)
    else:
        assert init['config']['layer'] == -1 and init['config']['head_width'] == A.HEAD_W


def new_head(cfg, seed=0):
    return MD.new_head(tuple(cfg['layers']), seed) if NT else MD.new_head(seed)


class Ctx:
    """Loaded encoder + what is needed to call it."""
    def __init__(self, enc, cfg, lora_api=None):
        self.enc, self.cfg, self.lora_api = enc, cfg, lora_api
        if NT:
            self.esm, self.vocab = enc
            self.layers, self.norm = tuple(cfg['layers']), cfg['layer_norm_per_layer']
        else:
            self.layer = cfg['layer']


def _pc2_lora_api():
    import importlib.util
    p = T.EXPS / 'M35-PLANTCAD2-LORA' / 'scripts' / 'pc2_lora.py'
    if not p.exists():
        raise RuntimeError(f'PlantCAD2 LoRA requested but {p} does not exist yet (M35 defines it); use --lora none for the frozen-encoder variant')
    spec = importlib.util.spec_from_file_location('pc2_lora', p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for f in ('add_lora', 'lora_state', 'load_lora_state', 'CFG'):
        assert hasattr(mod, f), f'pc2_lora.{f} missing'
    return mod


def load(cfg, lora=True, device='cuda'):
    """-> (Ctx, trainable LoRA params list)."""
    if NT:
        enc = A.load_ntv2(device)
        params = A.add_lora(enc[0], cfg['lora_encoder_layers'], cfg['lora_r'], cfg['lora_alpha']) if lora else []
        return Ctx(enc, cfg), params
    model = A.load_pc2(device)
    if not lora:
        return Ctx(model, cfg), []
    api = _pc2_lora_api()
    params = api.add_lora(model, api.CFG)
    cfg['lora'] = dict(api.CFG)
    return Ctx(model, cfg, api), params


def lora_state(ctx):
    return A.lora_state(ctx.esm) if NT else (ctx.lora_api.lora_state(ctx.enc) if ctx.lora_api else {})


def copy_lora(ctx, state):
    """Load a LoRA state into the ALREADY ADAPTED encoder (resume)."""
    if NT:
        for k, v in state.items():
            i, n, w = k.split('.')
            getattr(getattr(ctx.esm.encoder.layer[int(i)].attention.self, n), w).data.copy_(v)
    elif ctx.lora_api:
        ctx.lora_api.load_lora_state(ctx.enc, state, ctx.lora_api.CFG)


def apply_saved_lora(ctx, st):
    """Inference: adapt the encoder from a model file (`st` = torch.load of model_N.pt)."""
    if NT:
        if 'lora' in st:
            A.load_lora_state(ctx.esm, st['lora'], st['config']['lora_r'], st['config']['lora_alpha'])
            return True
        return False
    if st.get('lora'):
        api = _pc2_lora_api()
        ctx.lora_api = api
        api.add_lora(ctx.enc, st['config']['lora'])
        api.load_lora_state(ctx.enc, st['lora'], st['config']['lora'])
        return True
    return False


def load_for_inference(st, device='cuda'):
    cfg = st['config']
    enc = A.load_ntv2(device) if NT else A.load_pc2(device)
    ctx = Ctx(enc, cfg)
    apply_saved_lora(ctx, st)
    return ctx


def geometry(a, b):
    return A.nt_geometry(a, b) if NT else A.pc_geometry(a, b)


# ---- training items -------------------------------------------------------------------------------------------------------------------------------
def item_windows(chromosome, geo, strand):
    """-> dict(seqs=oriented window strings [, gen=genomic strings])"""
    if NT:
        return dict(seqs=MD.oriented_windows(chromosome, geo, strand))
    gen = MD.genomic_window_strings(chromosome, geo)
    return dict(seqs=MD.orient_strings(gen, strand), gen=gen)


def _pc2_encode(enc, strings, layer=-1, batch=8, grad=False, device='cuda'):
    """As M32 `encode_windows`, but differentiable when `grad` (the M32 one is inference_mode)."""
    out = [None] * len(strings)
    by_len = {}
    for i, s in enumerate(strings):
        if not 0 < len(s) <= A.PC_WIDTH:
            raise ValueError('Nonempty window of at most PC_WIDTH bp required')
        by_len.setdefault(len(s), []).append(i)
    with (torch.enable_grad() if grad else torch.no_grad()):
        for n, idx in by_len.items():
            for j in range(0, len(idx), batch):
                part = idx[j:j + batch]
                ids = torch.stack([A.token_ids(strings[i], device) for i in part])
                with torch.autocast('cuda', dtype=torch.bfloat16):
                    h = enc(ids).last_hidden_state if layer == -1 else enc(ids, output_hidden_states=True).hidden_states[layer]
                if h.shape != (len(part), n, A.HIDDEN):
                    raise ValueError(f'unexpected feature shape {tuple(h.shape)}')
                if not bool(torch.isfinite(h).all()):
                    raise ValueError('Finite PlantCAD2 features required')
                for k, i in enumerate(part):
                    out[i] = h[k].to(torch.bfloat16)
    return out


def window_pots(head, ctx, it, idx, grad, device='cuda'):
    """Valid-length 15-role potentials of windows `idx` of a draw item (list of [len,15] tensors on the GPU); autograd through encoder (LoRA) and head iff `grad`."""
    seqs = [it['seqs'][i] for i in idx]
    if NT:
        feats = A.encode_windows(ctx.esm, ctx.vocab, seqs, ctx.layers, batch=len(idx), norm=ctx.norm, custom=True, grad=grad)
    else:
        gen = _pc2_encode(ctx.enc, [it['gen'][i] for i in idx], ctx.layer, batch=8, grad=grad)
        feats = A.oriented_tokens(gen, it['strand'])
    return MD.window_potentials(head, feats, seqs, device=device)


# ---- inference ---------------------------------------------------------------------------------------------------------------------------------------
def region_potentials(head, ctx, sequence, strands, batch=4, device='cuda'):
    """{strand: oriented-region 15-role potentials [L,15] (CPU FP32)}, number of encoder windows.  Nearest-centre ownership via the frozen stitcher.
    ntv2: the M31 `stitched_potentials` per strand (streams the windows).  pc2: the M32 scheme -- one forward per genomic window serves both strands."""
    if NT:
        import c1_infer as I                       # M31 (unchanged)
        out, n = {}, 0
        for st in strands:
            out[st], n = I.stitched_potentials(head, ctx.esm, ctx.vocab, ctx.layers, ctx.norm, sequence, st, batch=batch, custom=True)
        return out, n
    from src.m28.contig_emission_lattice import stitch_contig_emissions
    import manifest as MF
    L = len(sequence)
    locations = [(a, min(a + A.PC_WIDTH, L)) for a in MF.starts(L, A.PC_WIDTH, A.PC_STRIDE)]
    recs = {st: [] for st in strands}
    with torch.inference_mode():
        for i in range(0, len(locations), batch):
            gen = [sequence[a:b] for a, b in locations[i:i + batch]]
            feats = A.encode_windows(ctx.enc, gen, ctx.layer, batch=batch)
            for st in strands:
                seqs = MD.orient_strings(gen, st)
                for j, (s, lg) in enumerate(zip(seqs, MD.window_potentials(head, A.oriented_tokens(feats, st), seqs, device=device))):
                    recs[st].append(dict(window_index=i + j, sequence=s, local_potentials=lg.float().cpu()))
            del feats
        out = {st: stitch_contig_emissions(sequence, locations, st, iter(recs[st]))['local_potentials'] for st in strands}
    return out, len(locations)
