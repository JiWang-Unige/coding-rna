"""M30 inference (copy of a1_infer with M30 contigs/offsets and the bias-free A2 head).
Original: M29-A1 inference: joint-24,576 frozen-encoder head logits + fixed Carbon-A offset -> exact skyline decode (frozen decoder, M29 arena copy).
GPU encodes contig/strand k+1 while CPU processes decode contig/strand k.  Resumable per (contig, strand): head logits are kept under
<out>/head/<seq>_<strand>.npy and chains under <out>/chains/<seq>_<strand>.json.
Usage: a1_infer.py --model model_N.pt --seqs NC_..,NC_.. --out DIR [--decode-workers 6] [--cpus 2,4,..]
       a1_infer.py --gate --seqs ...     (head = 0: the parity-gate model through THIS code path; sanity check / baseline)
"""
import argparse, json, multiprocessing as mp, os, re, sys, time
import numpy as np
import torch
import m30_common as M30
import a1_common as K
import a1_offset as O
M30.register()
O.offset_path = M30.offset_path      # M30 offsets (module-level so spawned decode workers see it too)


def _decode_init(cpus):
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(1)
    import a1_common  # noqa: F401
    a1_common.repoint()


def _decode_job(seq, strand, head_path, oriented_len, out_json, gate, bias_path=None, head_scale=1.0):
    import numpy as np
    from m29_skyline_reuse import decode_full_grammar_skyline
    from src.m28.c0_decode import genomic_chain
    tracks = O.load_tracks(seq)
    c = tracks[strand]
    L = c.shape[1]
    if gate:
        lg = np.zeros((L, 15), dtype=np.float32)
    else:
        lg = np.load(head_path)
        assert lg.shape == (L, 15) and lg.dtype == np.float32
        if head_scale != 1.0:
            lg *= np.float32(head_scale)      # diagnostic: shrink the residual head before adding the offset
        if bias_path:
            lg -= np.load(bias_path).astype(np.float32)[None, :]      # post-hoc diagnostic: remove a constant per-role head offset (estimated on TRAIN)
    step = 2_000_000
    for a in range(0, L, step):
        b = min(L, a + step)
        lg[a:b] += O.expand(c, a, b)
    seqstr = _read(seq)
    oriented = seqstr if strand == '+' else O.rc(seqstr)
    t0 = time.perf_counter()
    res = decode_full_grammar_skyline(oriented, torch.from_numpy(lg), 128)
    sec = time.perf_counter() - t0
    chains = []
    for chain, gain in zip(res['chains'], res['chain_gains']):
        g = genomic_chain(chain, 0, L, strand)
        chains.append(dict(exons=[(int(a), int(b)) for a, b in g], gain=float(gain)))
    json.dump(dict(seq=seq, strand=strand, chains=chains, total_gain=res['gain'], decode_seconds=sec), open(out_json, 'w'))
    return out_json


def _read(seq):
    import m29_common as M
    return re.sub('[^ACGT]', 'N', M.read_sequence(seq))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model')
    ap.add_argument('--gate', action='store_true')
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--decode-workers', type=int, default=6)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--keep-head', action='store_true', default=True)
    ap.add_argument('--head-from', default=None, help='reuse head logits of another run dir (no GPU)')
    ap.add_argument('--bias', default=None)
    ap.add_argument('--head-scale', type=float, default=1.0)
    args = ap.parse_args()
    out = K.Path(args.out) if hasattr(K, 'Path') else __import__('pathlib').Path(args.out)
    for sub in ('head', 'chains'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    cpus = None
    if args.cpus:
        cpus = []
        for part in args.cpus.split(','):
            lo, _, hi = part.partition('-')
            cpus += list(range(int(lo), int(hi or lo) + 1))
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    seqs = args.seqs.split(',')
    pool = mp.get_context('spawn').Pool(args.decode_workers, initializer=_decode_init, initargs=(cpus,))
    jobs = []
    began = time.perf_counter()
    if not args.gate and not args.head_from:
        import a1_model as MD
        import a2_model as MD2
        import annotate_full_grammar_joint24 as A
        import annotate_full_grammar_coarse_channel as R
        K.repoint()
        st = torch.load(args.model, map_location='cpu', weights_only=False)
        model = MD2.new_model(0)
        model.load_state_dict(st['model'])
        model = model.eval().requires_grad_(False).to('cuda')
        backbone, tokenizer = MD.G.load_native()
        torch.cuda.reset_peak_memory_stats()
    enc_s = 0.
    for seq in seqs:
        seqstr = _read(seq)
        L = len(seqstr)
        locations = None
        for strand in '+-':
            cj = out / 'chains' / f'{seq}_{strand}.json'
            hp = (K.Path(args.head_from) if args.head_from else out) / 'head' / f'{seq}_{strand}.npy'
            if cj.exists():
                continue
            if not args.gate and not hp.exists() and not args.head_from:
                import annotate_full_grammar_joint24 as A
                import annotate_full_grammar_coarse_channel as R
                locations = locations or [(a, min(a + R.M.WIDTH, L)) for a in R.M.starts(L)]
                t0 = time.perf_counter()
                joined = A.contig_emissions(model, backbone, tokenizer, seqstr, locations, strand, 'cuda')
                np.save(str(hp) + '.tmp.npy', joined['local_potentials'].numpy())
                os.replace(str(hp) + '.tmp.npy', hp)
                enc_s += time.perf_counter() - t0
                del joined
            jobs.append(pool.apply_async(_decode_job, (seq, strand, str(hp), L, str(cj), args.gate, args.bias, args.head_scale)))
            print(json.dumps(dict(event='submitted_decode', seq=seq, strand=strand, t=time.perf_counter() - began)), flush=True)
    for j in jobs:
        print(json.dumps(dict(event='decoded', f=j.get(), t=time.perf_counter() - began)), flush=True)
    pool.close(); pool.join()
    # GFF3 per contig
    from src.m28.export import as_transcript, write_gff3
    import m29_common as M
    for seq in seqs:
        L = M.LENGTH[seq]
        txs = []
        for strand in '+-':
            d = json.load(open(out / 'chains' / f'{seq}_{strand}.json'))
            for i, c in enumerate(d['chains']):
                txs.append(as_transcript([tuple(e) for e in c['exons']], seq, strand, f'M29_A1_{seq}_{"p" if strand=="+" else "m"}_{i}', c['gain']))
        (out / seq / 'prediction').mkdir(parents=True, exist_ok=True)
        gff = out / seq / 'prediction' / 'full_grammar.gff3'
        if gff.exists():
            gff.unlink()
        write_gff3(gff, txs, {seq: L}, 'M29_A1')
        print(json.dumps(dict(event='gff', seq=seq, chains=len(txs))), flush=True)
    json.dump(dict(completed=True, seqs=seqs, model=args.model, gate=args.gate, encode_seconds=enc_s, wall_seconds=time.perf_counter() - began,
                   peak_gpu_gb=(torch.cuda.max_memory_allocated() / 2**30 if not (args.gate or args.head_from) else None)), open(out / f'DONE_{"_".join(seqs)[:60]}.json', 'w'))


if __name__ == '__main__':
    main()
