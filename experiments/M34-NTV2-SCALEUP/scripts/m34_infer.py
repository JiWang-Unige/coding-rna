"""M34 inference = M31 `c1_infer.py` (NT-v2 windows + LoRA + C0 head -> whole-region exact skyline decode -> GFF3) with the new-species validation
pseudo-contigs registered in the frozen M29 helpers (also inside the spawned decode workers).  The M31 file is imported, not edited; the model
loading, stitching and decoding code is the same.  Regions are SEQ or SEQ:START-END, GFF3 coordinates are global.
Usage: m34_infer.py --model model_N.pt --regions A,B,... --out DIR [--decode-workers 4] [--cpus ...]
"""
import argparse, hashlib, json, multiprocessing as mp, os, re, time
import numpy as np
import torch

import m34_train_common as T
import m34_backend as BK
import c1_infer as I                    # M31 or M32 (unchanged; only the encoder-independent helpers parse_region / resolve / decode are used)
M = T.M


def _decode_init(cpus):
    I._decode_init(cpus)
    T.register()


def _decode_job(seq, a, b, strand, pot_path, out_json, keep):
    T.register()
    return I._decode_job(seq, a, b, strand, pot_path, out_json, keep)


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
    MC = T.register()
    regions = [I.parse_region(r) for r in args.regions.split(',')]
    st = torch.load(args.model, map_location='cpu', weights_only=False)
    cfg = st['config']
    model = BK.new_head(cfg)
    model.load_state_dict(st['model'])
    model = model.eval().requires_grad_(False).to('cuda')
    ctx = BK.load_for_inference(st)
    torch.cuda.reset_peak_memory_stats()
    pool = mp.get_context('spawn').Pool(args.decode_workers, initializer=_decode_init, initargs=(cpus,))
    jobs, began, enc_s, bases, nwin = [], time.perf_counter(), 0., 0, 0
    for seq, a, b in regions:
        full = re.sub('[^ACGT]', 'N', MC.read_sequence(seq))
        b, name = I.resolve(seq, a, b, len(full))
        sl = full[a:b]
        todo = [st_ for st_ in '+-' if not (out / 'chains' / f'{name}_{st_}.json').exists()]
        if not todo:
            continue
        t0 = time.perf_counter()
        pots, n = BK.region_potentials(model, ctx, sl, todo)            # ntv2: per strand (M31 code); pc2: one forward serves both strands
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
        _, name = I.resolve(seq, a, b, L)
        txs = []
        for strand in '+-':
            d = json.load(open(out / 'chains' / f'{name}_{strand}.json'))
            for i, c in enumerate(d['chains']):
                txs.append(as_transcript([tuple(e) for e in c['exons']], seq, strand, f'M34{T.BACKEND}_{name}_{"p" if strand == "+" else "m"}_{i}', c['gain']))
        (out / name / 'prediction').mkdir(parents=True, exist_ok=True)
        gff = out / name / 'prediction' / 'full_grammar.gff3'
        if gff.exists():
            gff.unlink()
        write_gff3(gff, txs, {seq: L}, 'M34_NTV2')
        print(json.dumps(dict(event='gff', region=name, chains=len(txs))), flush=True)
    wall = time.perf_counter() - began
    json.dump(dict(completed=True, regions=args.regions, model=args.model, step=st.get('step'), backend=T.BACKEND, encode_seconds=enc_s, wall_seconds=wall,
                   bases_per_strand=bases, encoder_windows=nwin, encode_kbp_per_s_both_strands=2 * bases / 1000 / max(enc_s, 1e-9),
                   wall_kbp_per_s_both_strands=2 * bases / 1000 / wall, peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30,
                   peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30), open(out / f'DONE_{hashlib.md5(args.regions.encode()).hexdigest()[:8]}.json', 'w'))


if __name__ == '__main__':
    main()
