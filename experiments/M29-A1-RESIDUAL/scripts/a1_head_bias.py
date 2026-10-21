"""Mean head logit per role over a TRAIN sample (post-hoc diagnostic vector for --bias)."""
import sys, numpy as np, os
O = os.path.expanduser('~/coding-rna/outputs/M29-A1-RESIDUAL')
run = sys.argv[1]
acc, n = np.zeros(15), 0
for seq in ('NC_003070.9', 'NC_003071.7', 'NC_003076.8', 'NC_089043.1', 'NC_089044.1', 'NC_089046.1'):
    for st in '+-':
        h = np.load(f'{run}/head/{seq}_{st}.npy', mmap_mode='r')
        x = np.asarray(h[2_000_000:8_000_000], dtype=np.float64)
        acc += x.sum(0); n += len(x)
np.save(sys.argv[2], acc / n)
print(np.round(acc / n, 3))
