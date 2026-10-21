"""Parse every tool's TRAIN GFF once into compact per-contig records (cached pickle)."""
import pickle
import sys
import m29_common as M
from m29_r3_posterior_filter import carbon_conf


def tool_chains(tool, seq):
    p = M.CACHE / f'pred_{tool}_{seq}.pkl'
    if p.exists():
        return pickle.loads(p.read_bytes())
    gff = M.gff_path(tool, seq)
    conf = carbon_conf(gff) if tool.startswith('CarbonA') else None
    out = {}
    for tx in M.load_txs(gff, seq):
        k = M.key(tx)
        if k in out:
            continue
        out[k] = dict(key=k, strand=tx['strand'], exons=[(a, b) for a, b, _ in tx['CDS']], id=tx['id'],
                      conf=(conf[tx['gene_id']] if conf is not None else 0.0), na=(conf is not None and conf[tx['gene_id']] is None))
    p.write_bytes(pickle.dumps(out, protocol=4))
    return out


if __name__ == '__main__':
    tools = sys.argv[1].split(',')
    for seq in sys.argv[2:]:
        for t in tools:
            c = tool_chains(t, seq)
            print(t, seq, len(c), flush=True)
