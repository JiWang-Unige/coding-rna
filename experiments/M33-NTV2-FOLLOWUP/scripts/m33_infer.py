"""M33-NTV2-FOLLOWUP inference: copy of M31 `c1_infer.py` (frozen M31 file untouched) with three changes only:
  (1) M30 contigs (Dryas / Eschscholzia) are registered in the M29 helpers, also inside the spawned decode workers;
  (2) the oriented 15-role potentials  pot/<name>_<strand>.npy  are KEPT (needed for the gene-set posterior, `m33_posterior.py`);
  (3) the output GFF3 id prefix is M33.
Model B (step 1536: C0 head + NT-v2 q/v LoRA, hidden_states 21+27), 12,288-bp windows / 6,144 stride, exact whole-contig skyline decode.
Usage: m33_infer.py --model model_001536.pt --regions SEQ,SEQ:START-END,... --out DIR [--decode-workers 4] [--cpus 2,4,...]
"""
import argparse, hashlib, json, multiprocessing as mp, os, re, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[0].parent                               # .../experiments
for p in (str(EXP / 'M31-C1-NTV2' / 'scripts'), str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)
import numpy as np
import torch

import ntv2_adapter as A
import m31_common as M
import c1_model as MD
import manifest as MF
import m30_common as M30
M30.register()                                             # module level: also runs in spawned workers (__mp_main__)


def parse_region(r):
    seq, _, span = r.partition(':')
    if not span:
        return seq, 0, None
    a, _, b = span.partition('-')
    return seq, int(a), int(b)


def resolve(seq, a, b, L):
    b = L if b is None else min(b, L)
    return b, (seq if (a == 0 and b == L) else f'{seq}_{a}_{b}')


def _decode_init(cpus):
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(1)
    import a1_common
    a1_common.repoint()
    M30.register()


def _decode_job(seq, a, b, strand, pot_path, out_json, keep):
    from m29_skyline_reuse import decode_full_grammar_skyline
    from src.m28.c0_decode import genomic_chain
    import m29_common as MC
    lg = np.load(pot_path)
    full = re.sub('[^ACGT]', 'N', MC.read_sequence(seq))
    sl = full[a:b]
    L = len(sl)
    assert lg.shape == (L, 15) and lg.dtype == np.float32
    oriented = sl if strand == '+' else MF.rc(sl)
    t0 = time.perf_counter()
    res = decode_full_grammar_skyline(oriented, torch.from_numpy(lg), 128)
    sec = time.perf_counter() - t0
    chains = []
    for chain, gain in zip(res['chains'], res['chain_gains']):
        g = genomic_chain(chain, 0, L, strand)
        chains.append(dict(exons=[(int(x) + a, int(y) + a) for x, y in g], gain=float(gain)))
    json.dump(dict(seq=seq, start=a, end=b, strand=strand, chains=chains, total_gain=res['gain'], decode_seconds=sec), open(out_json, 'w'))
    if not keep:
        os.remove(pot_path)
    return out_json


def stitched_potentials(model, esm, vocab, layers, norm, sequence, strand, batch=4, custom=False):
    """Oriented-region 15-role potentials [L,15] (CPU FP32), nearest-centre owned, via the frozen stitcher (identical to M31)."""
    from src.m28.contig_emission_lattice import stitch_contig_emissions
    L = len(sequence)
    locations = [(a, min(a + A.NT_WIDTH, L)) for a in MF.starts(L, A.NT_WIDTH, A.NT_STRIDE)]

    def records():
        for i in range(0, len(locations), batch):
            chunk = locations[i:i + batch]
            seqs = [sequence[a:b] if strand == '+' else MF.rc(sequence[a:b]) for a, b in chunk]
            feats = A.encode_windows(esm, vocab, seqs, layers, batch=batch, norm=norm, custom=custom)
            pots = MD.window_potentials(model, feats, seqs)
            for j, (s, lg) in enumerate(zip(seqs, pots)):
                yield dict(window_index=i + j, sequence=s, local_potentials=lg.float().cpu())
    with torch.inference_mode():
        return stitch_contig_emissions(sequence, locations, strand, records())['local_potentials'], len(locations)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--regions', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--decode-workers', type=int, default=4)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--drop-pot', action='store_true', help='delete potentials after decoding (default: keep, needed for posteriors)')
    args = ap.parse_args()
    out = Path(args.out)
    for sub in ('pot', 'chains'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    import m29_common as MC
    regions = [parse_region(r) for r in args.regions.split(',')]
    st = torch.load(args.model, map_location='cpu', weights_only=False)
    cfg = st['config']
    layers, norm = tuple(cfg['layers']), cfg['layer_norm_per_layer']
    model = MD.new_head(layers)
    model.load_state_dict(st['model'])
    model = model.eval().requires_grad_(False).to('cuda')
    esm, vocab = A.load_ntv2()
    custom = 'lora' in st
    if custom:
        A.load_lora_state(esm, st['lora'], cfg['lora_r'], cfg['lora_alpha'])
    torch.cuda.reset_peak_memory_stats()
    pool = mp.get_context('spawn').Pool(args.decode_workers, initializer=_decode_init, initargs=(cpus,))
    jobs, began, enc_s, bases, nwin = [], time.perf_counter(), 0., 0, 0
    for seq, a, b in regions:
        full = re.sub('[^ACGT]', 'N', MC.read_sequence(seq))
        b, name = resolve(seq, a, b, len(full))
        sl = full[a:b]
        for strand in '+-':
            cj = out / 'chains' / f'{name}_{strand}.json'
            if cj.exists():
                continue
            t0 = time.perf_counter()
            pp = out / 'pot' / f'{name}_{strand}.npy'
            if pp.exists():                                   # resume: potentials already encoded, only decode again
                pot = None
            else:
                pot, n = stitched_potentials(model, esm, vocab, layers, norm, sl, strand, custom=custom)
                np.save(str(pp) + '.tmp.npy', pot.numpy())
                os.replace(str(pp) + '.tmp.npy', pp)
                enc_s += time.perf_counter() - t0
                bases += len(sl); nwin += n
            del pot
            jobs.append(pool.apply_async(_decode_job, (seq, a, b, strand, str(pp), str(cj), not args.drop_pot)))
            print(json.dumps(dict(event='submitted_decode', region=name, strand=strand, encode_s=round(time.perf_counter() - t0, 1), t=round(time.perf_counter() - began, 1))), flush=True)
    for j in jobs:
        print(json.dumps(dict(event='decoded', f=j.get(), t=round(time.perf_counter() - began, 1))), flush=True)
    pool.close(); pool.join()
    from src.m28.export import as_transcript, write_gff3
    for seq, a, b in regions:
        L = MC.LENGTH[seq]
        _, name = resolve(seq, a, b, L)
        txs = []
        for strand in '+-':
            d = json.load(open(out / 'chains' / f'{name}_{strand}.json'))
            for i, c in enumerate(d['chains']):
                txs.append(as_transcript([tuple(e) for e in c['exons']], seq, strand, f'M33_{name}_{"p" if strand == "+" else "m"}_{i}', c['gain']))
        (out / name / 'prediction').mkdir(parents=True, exist_ok=True)
        gff = out / name / 'prediction' / 'full_grammar.gff3'
        if gff.exists():
            gff.unlink()
        write_gff3(gff, txs, {seq: L}, 'M33_NTV2_B')
        print(json.dumps(dict(event='gff', region=name, chains=len(txs))), flush=True)
    wall = time.perf_counter() - began
    json.dump(dict(completed=True, regions=args.regions, model=args.model, step=st.get('step'), layers=list(layers), encode_seconds=enc_s, wall_seconds=wall,
                   bases_per_strand=bases, nt_windows=nwin, encode_kbp_per_s_both_strands=2 * bases / 1000 / max(enc_s, 1e-9),
                   wall_kbp_per_s_both_strands=2 * bases / 1000 / wall, peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30,
                   peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30), open(out / f'DONE_{hashlib.md5(args.regions.encode()).hexdigest()[:8]}.json', 'w'))


if __name__ == '__main__':
    main()
