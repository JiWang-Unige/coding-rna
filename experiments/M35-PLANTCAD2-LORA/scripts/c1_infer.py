"""M35-PLANTCAD2-LORA inference: PlantCAD2-Small (frozen weights + the LoRA stored in the checkpoint, if any) windows (8,192 bp, stride 6,144, both strands from one forward, nearest-centre ownership via the frozen stitcher) + trained C0 head
-> whole-region exact skyline decode (frozen decoder, M29 arena copy) -> GFF3.  No offset, no external model, no reference.
GPU encodes (region, strand) k+1 while CPU workers decode k.  Resumable per (region, strand): chains/<name>_<strand>.json.
Region = (seq, start, end) in genomic coordinates, spelled  SEQ  or  SEQ:START-END ; GFF3 coordinates are global.
Usage: c1_infer.py --model model_N.pt --regions NC_003071.7:5000000-13000000,NC_089044.1:14000000-26000000 --out DIR [--decode-workers 4]
"""
import argparse, hashlib, json, multiprocessing as mp, os, re, sys, time
import numpy as np
import torch

import pc2_adapter as A
import pc2_lora as PL
import m35_common as M
import c1_model as MD
import manifest as MF


def parse_region(r):
    seq, _, span = r.partition(':')
    if not span:
        return seq, 0, None
    a, _, b = span.partition('-')
    return seq, int(a), int(b)


def resolve(seq, a, b, L):
    """-> (end clamped to the contig, file-name stem)"""
    b = L if b is None else min(b, L)
    return b, (seq if (a == 0 and b == L) else f'{seq}_{a}_{b}')


def _decode_init(cpus):
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(1)
    import a1_common
    a1_common.repoint()


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


def region_potentials(model, enc, layer, sequence, batch=4, device='cuda'):
    """Both strands' oriented-region 15-role potentials ([L,15] CPU FP32 each), nearest-centre owned via the frozen stitcher.
    PlantCAD2 is exactly RC-equivariant, so ONE forward over each genomic window gives the features of both strands."""
    from src.m28.contig_emission_lattice import stitch_contig_emissions
    L = len(sequence)
    locations = [(a, min(a + A.PC_WIDTH, L)) for a in MF.starts(L, A.PC_WIDTH, A.PC_STRIDE)]
    recs = {'+': [], '-': []}
    with torch.inference_mode():
        for i in range(0, len(locations), batch):
            chunk = locations[i:i + batch]
            gen = [sequence[a:b] for a, b in chunk]
            feats = A.encode_windows(enc, gen, layer, batch=batch, device=device)
            for strand in '+-':
                seqs = MD.orient_strings(gen, strand)
                pots = MD.window_potentials(model, A.oriented_tokens(feats, strand), seqs, device=device)
                for j, (s, lg) in enumerate(zip(seqs, pots)):
                    recs[strand].append(dict(window_index=i + j, sequence=s, local_potentials=lg.float().cpu()))
            del feats
        out = {st: stitch_contig_emissions(sequence, locations, st, iter(recs[st]))['local_potentials'] for st in '+-'}
    return out, len(locations)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True)
    ap.add_argument('--regions', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--decode-workers', type=int, default=4)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--keep-head', action='store_true')
    args = ap.parse_args()
    out = M.Path(args.out)
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
    layer = cfg['layer']
    model = MD.new_head()
    model.load_state_dict(st['model'])
    model = model.eval().requires_grad_(False).to('cuda')
    enc = A.load_pc2()
    if 'lora' in st:                      # M35: frozen PlantCAD2 + the LoRA of the checkpoint (identity-at-init adapters; patched Mamba2 forward, no grad)
        PL.load_lora_state(enc, st['lora'], cfg.get('lora_r', PL.RANK), cfg.get('lora_alpha', PL.ALPHA))
    torch.cuda.reset_peak_memory_stats()
    pool = mp.get_context('spawn').Pool(args.decode_workers, initializer=_decode_init, initargs=(cpus,))
    jobs, began, enc_s, bases, nwin = [], time.perf_counter(), 0., 0, 0
    for seq, a, b in regions:
        full = re.sub('[^ACGT]', 'N', MC.read_sequence(seq))
        b, name = resolve(seq, a, b, len(full))
        sl = full[a:b]
        todo = [st for st in '+-' if not (out / 'chains' / f'{name}_{st}.json').exists()]
        if not todo:
            continue
        t0 = time.perf_counter()
        pots, n = region_potentials(model, enc, layer, sl)
        enc_s += time.perf_counter() - t0
        bases += len(sl); nwin += n
        for strand in todo:
            cj = out / 'chains' / f'{name}_{strand}.json'
            pp = out / 'pot' / f'{name}_{strand}.npy'
            np.save(str(pp) + '.tmp.npy', pots[strand].numpy())
            os.replace(str(pp) + '.tmp.npy', pp)
            jobs.append(pool.apply_async(_decode_job, (seq, a, b, strand, str(pp), str(cj), args.keep_head)))
            print(json.dumps(dict(event='submitted_decode', region=name, strand=strand, encode_s=round(time.perf_counter() - t0, 1), t=round(time.perf_counter() - began, 1))), flush=True)
        del pots
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
                txs.append(as_transcript([tuple(e) for e in c['exons']], seq, strand, f'M35_{name}_{"p" if strand == "+" else "m"}_{i}', c['gain']))
        (out / name / 'prediction').mkdir(parents=True, exist_ok=True)
        gff = out / name / 'prediction' / 'full_grammar.gff3'
        if gff.exists():
            gff.unlink()
        write_gff3(gff, txs, {seq: L}, 'M35_PLANTCAD2_LORA')
        print(json.dumps(dict(event='gff', region=name, chains=len(txs))), flush=True)
    wall = time.perf_counter() - began
    json.dump(dict(completed=True, regions=args.regions, model=args.model, step=st.get('step'), layer=layer, lora_blocks=sorted({int(k.split('.')[0]) for k in st.get('lora', {})}), encode_seconds=enc_s, wall_seconds=wall,
                   bases_per_strand=bases, pc2_windows=nwin, encode_kbp_per_s_both_strands=bases / 1000 / max(enc_s, 1e-9),
                   wall_kbp_per_s_both_strands=bases / 1000 / wall, peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30,
                   peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30), open(out / f'DONE_{hashlib.md5(args.regions.encode()).hexdigest()[:8]}.json', 'w'))


if __name__ == '__main__':
    main()
