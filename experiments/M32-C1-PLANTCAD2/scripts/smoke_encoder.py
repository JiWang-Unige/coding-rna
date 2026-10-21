"""PlantCAD2-Small encoder smoke (<= ~3 min GPU, ~3 GB): load, RC equivariance under autocast (flip-both vs swap-halves), B=1 vs B=4 consistency,
throughput / memory on real Arabidopsis TRAIN sequence, hidden_states count.  Writes JSON to stdout and --out.  Aborts unless >= 8 GB GPU memory is free."""
import argparse, json, subprocess, time
import numpy as np
import torch
import pc2_adapter as A
import m32_common as M
import manifest as MF


def free_mib():
    return int(subprocess.run(['nvidia-smi', '--query-gpu=memory.free', '--format=csv,noheader,nounits'], capture_output=True, text=True).stdout.split()[0])


def cos(a, b):
    c = torch.nn.functional.cosine_similarity(a.float(), b.float(), dim=-1)
    return float(c.mean()), float(c.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    if free_mib() < 8000:
        raise SystemExit(f'only {free_mib()} MiB free: not running (would risk the other GPU jobs)')
    import m29_common as MC
    seq = MC.read_sequence('NC_003070.9')[3_000_000:3_000_000 + 8192 * 8]
    wins = [seq[i * 8192:(i + 1) * 8192] for i in range(8)]
    res = dict(free_mib_before=free_mib())
    t0 = time.perf_counter()
    enc = A.load_pc2()
    res['load_s'] = round(time.perf_counter() - t0, 1)
    res['n_params'] = sum(p.numel() for p in enc.parameters())
    res['weights_dtype'] = str({p.dtype for p in enc.parameters()})
    torch.cuda.reset_peak_memory_stats()
    f = A.encode_windows(enc, wins[:2], batch=2)
    fr = A.encode_windows(enc, [MF.rc(w) for w in wins[:2]], batch=2)
    d = A.D_MODEL
    r = {}
    for k in range(2):
        a = fr[k]
        both = f[k].flip(0).flip(-1)
        swap = torch.cat([f[k][..., d:], f[k][..., :d]], -1).flip(0)
        r[k] = dict(cos_flipboth=cos(a, both), cos_swaphalves=cos(a, swap), maxabs_flipboth=float((a.float() - both.float()).abs().max()))
    res['rc_equivariance_autocast'] = r
    single = A.encode_windows(enc, wins[:4], batch=1)
    batch4 = A.encode_windows(enc, wins[:4], batch=4)
    res['b1_vs_b4_cos'] = [cos(x, y) for x, y in zip(single, batch4)]
    res['feature_stats'] = dict(std=float(f[0].float().std()), absmax=float(f[0].float().abs().max()))
    with torch.inference_mode(), torch.autocast('cuda', dtype=torch.bfloat16):
        res['n_hidden_states'] = len(enc(A.token_ids(wins[0])[None], output_hidden_states=True).hidden_states)
    ts = []
    for _ in range(3):
        torch.cuda.synchronize(); t = time.perf_counter()
        A.encode_windows(enc, wins, batch=4)
        torch.cuda.synchronize(); ts.append(time.perf_counter() - t)
    res['encode_8x8192_B4_s'] = [round(x, 3) for x in ts]
    res['encode_kbp_per_s'] = round(8 * 8.192 / min(ts), 1)
    res['peak_alloc_GB'] = round(torch.cuda.max_memory_allocated() / 2**30, 2)
    res['free_mib_after'] = free_mib()
    s = json.dumps(res, indent=1)
    print(s)
    if args.out:
        open(args.out, 'w').write(s)


if __name__ == '__main__':
    main()
