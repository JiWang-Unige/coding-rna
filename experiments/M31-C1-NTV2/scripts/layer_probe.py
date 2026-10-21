"""Cheap layer probe (decides which frozen NT-v2 layer(s) feed the C0 head).  One linear probe per hidden-state index (0..29) on
the per-token 3-class region label (background / CDS / intron, label of the token's centre base, known positions only),
trained on windows of the four TRAIN chromosomes, scored on windows of the two internal-validation chromosomes.  No DEV.
Usage: layer_probe.py [--train 80] [--val 80] [--epochs 3]  -> OUT/layer_probe.json
"""
import argparse, json, time
import numpy as np
import torch
import m31_common as M
import ntv2_adapter as A


def pick(windows, rng, n, chroms, per_species=None):
    pool = [w for w in windows.values() if w['seqid'] in chroms and w['valid_bases'] == 24576]
    rng.shuffle(pool)
    out, cnt = [], {}
    for w in pool:
        if cnt.get(w['species'], 0) >= n // 2:
            continue
        out.append(w); cnt[w['species']] = cnt.get(w['species'], 0) + 1
        if len(out) >= n:
            break
    return out


def tokens(w, sequences, catalogs, esm, vocab):
    """-> (feats list over layers [2 windows concatenated: 4096 x 1024 each], labels [4096], keep mask)"""
    chromosome, geometry, context, offset, targets, counts = M.example(w, sequences, catalogs)
    core = context[offset:offset + 24576]
    region = np.asarray(targets['region'])[:24576]
    cen = region[3::6]                                   # centre base of every 6-mer
    halves = [core[:A.NT_WIDTH], core[A.NT_WIDTH:]]
    ids = torch.tensor([vocab.ids(s) for s in halves], device='cuda')
    with torch.inference_mode(), torch.autocast('cuda', dtype=torch.bfloat16):
        hs = esm(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True).hidden_states
    feats = [h[:, 1:].reshape(-1, A.HIDDEN) for h in hs]
    lab = torch.as_tensor(cen, dtype=torch.long, device='cuda')
    return feats, lab, (lab >= 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train', type=int, default=80)
    ap.add_argument('--val', type=int, default=80)
    ap.add_argument('--epochs', type=int, default=3)
    args = ap.parse_args()
    t0 = time.time()
    windows, strata, catalogs, original = M.PP.load_population()
    sequences, _ = M.C.P.read_train_DNA()
    rng = np.random.default_rng(0)
    tr = pick(windows, rng, args.train, set(M.TRAIN_SEQS) - set(M.HELDOUT))
    va = pick(windows, rng, args.val, set(M.HELDOUT))
    esm, vocab = A.load_ntv2()
    nl = A.N_LAYERS + 1
    heads = torch.nn.ModuleList([torch.nn.Linear(A.HIDDEN, 3) for _ in range(nl)]).cuda()
    opt = torch.optim.AdamW(heads.parameters(), lr=2e-3, weight_decay=1e-4)
    ln = lambda x: torch.nn.functional.layer_norm(x.float(), (A.HIDDEN,))
    print('setup_s', round(time.time() - t0, 1), flush=True)
    for ep in range(args.epochs):
        order = rng.permutation(len(tr))
        tot = 0.
        for i in order:
            feats, lab, keep = tokens(tr[i], sequences, catalogs, esm, vocab)
            if keep.sum() == 0:
                continue
            loss = sum(torch.nn.functional.cross_entropy(heads[k](ln(f[keep])), lab[keep]) for k, f in enumerate(feats))
            opt.zero_grad(); loss.backward(); opt.step(); tot += float(loss) / nl
        print(json.dumps(dict(epoch=ep + 1, mean_loss=round(tot / len(tr), 4), t=round(time.time() - t0))), flush=True)
    conf = torch.zeros(nl, 3, 3, dtype=torch.long)
    n_tok = 0
    with torch.no_grad():
        for w in va:
            feats, lab, keep = tokens(w, sequences, catalogs, esm, vocab)
            n_tok += int(keep.sum())
            for k, f in enumerate(feats):
                p = heads[k](ln(f[keep])).argmax(-1)
                conf[k] += torch.bincount(lab[keep] * 3 + p, minlength=9).reshape(3, 3).cpu()
    res = []
    for k in range(nl):
        c = conf[k].double()
        f1 = [float(2 * c[j, j] / (c[j].sum() + c[:, j].sum()).clamp(min=1)) for j in range(3)]
        res.append(dict(hidden_state=k, acc=float(c.diag().sum() / c.sum()), macro_F1=float(np.mean(f1)), F1_bg_cds_intron=[round(x, 4) for x in f1]))
    out = dict(train_windows=len(tr), val_windows=len(va), val_tokens=n_tok, epochs=args.epochs, layers=res,
               best_by_macro_F1=sorted(res, key=lambda r: -r['macro_F1'])[:6])
    M.OUT.mkdir(parents=True, exist_ok=True)
    (M.OUT / 'layer_probe.json').write_text(json.dumps(out, indent=1))
    for r in res:
        print(r['hidden_state'], round(r['acc'], 4), round(r['macro_F1'], 4), r['F1_bg_cds_intron'])
    print('seconds', round(time.time() - t0))


if __name__ == '__main__':
    main()
