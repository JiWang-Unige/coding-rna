"""Parity gate: head = 0, fixed offset built from Carbon-A's native probabilities / transitions, decoded by OUR
full-grammar exact skyline (src.m28.full_grammar_skyline).  Expected: reproduce Carbon-A UNFILTERED (TRAIN pooled ~76.06).

Offset construction (oriented strand, per base i; role index as in src/m28/full_grammar_boundary_messages.py):
  gain[r] = row[r] - row[0]  (row[0] = background = 0).  Carbon-A path score =
      sum_exon-like bases log p + sum_other bases log(1-p) + transition log-odds + intron prior,
  so with  d_i = log p_i - log(1-p_i)  (p clipped to [1e-6, 1-1e-6], float16 input as in Carbon-A):
    roles with coding bases (1,2,3,7..14) : d_i          roles 4,5,6 (intron) and 0 (background): 0
    enter_gene  (first base of ATG)             : +tr[i,0]            on role 13 at an 'A' base starting ATG
    leave_gene  (last base of stop codon)       : +tr[j,1]            on role 14 at j when seq[j-2:j+1] in {TAA,TAG,TGA}
    enter_intron(first intron base e+1)         : +tr[e+1,2]+prior    on roles 7,8,9 (last exon base e) when seq[e+1:e+3] in {GT,GC}
    leave_intron(first exon base after intron)  : +tr[i,3]            on roles 10,11,12
    GC-AG prior log(0.009/0.99); GT-AG 0; AT-AC not representable in our grammar.
Variants: A_full (emission+transitions+prior), B_emission (emission only), C_trans_only_gene? (not run)
Region mode (smoke): decode only [a,b) of a contig and compare with Carbon-A post-processing of the SAME tracks.
"""
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CA_SRC = HERE.parent / 'carbon_a_pipeline_copy' / 'src'
sys.path.insert(0, str(CA_SRC))
import torch  # noqa: E402
NT = int(os.environ.get('M29_THREADS', '5'))
torch.set_num_threads(NT)
import m29_common as M  # noqa: E402  (also puts pipeline/m28 on sys.path)
from constants import (CODON_TABLE_PATH, TRANSITION_NETWORK_PATH, CONFIDENCE_NETWORK_PATH, MIN_PROB, MAX_PROB,  # noqa: E402
                       GC_AG_INTRON_LOG_PRIOR, POSITIVE_STRAND_PROBABILITY_COLUMN, NEGATIVE_STRAND_PROBABILITY_COLUMN,
                       RECORD_ID_COLUMN, SEQUENCE_COLUMN, MAX_BOUNDARY_INTRON_LENGTH_BP, MIN_CDS_LENGTH_BP)
from utils import probability_vector_from_parquet, encode_both_strands, load_codon_table, stop_codons_from_table  # noqa: E402
from transition_network import load_transition_network  # noqa: E402
from src.m28.full_grammar_skyline import decode_full_grammar_skyline  # noqa: E402
from src.m28.c0_decode import genomic_chain  # noqa: E402
from src.m28.export import as_transcript, write_gff3  # noqa: E402

OUTROOT = M.M29 / 'R6'
CACHE = M.CACHE
RAW = M.OUTS / 'M28-CARBON-A-TRAIN-RAW'
STOPS = ('TAA', 'TAG', 'TGA')
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def load_tracks(seq):
    import pyarrow.parquet as pq
    p = sorted((RAW / seq / 'native').glob('*.raw_probs.parquet'))[0]
    pf = pq.ParquetFile(p)
    batches = list(pf.iter_batches(batch_size=1, columns=[RECORD_ID_COLUMN, SEQUENCE_COLUMN, POSITIVE_STRAND_PROBABILITY_COLUMN,
                                                           NEGATIVE_STRAND_PROBABILITY_COLUMN]))
    assert len(batches) == 1
    b = batches[0]
    sequence = b.column(1)[0].as_py().upper()
    pos = probability_vector_from_parquet(b.column(2)[0])[:len(sequence)]
    neg = probability_vector_from_parquet(b.column(3)[0])[:len(sequence)]
    return b.column(0)[0].as_py(), sequence, pos, neg


def transition_tracks(tn, probs, key):
    f = CACHE / f'tt_{key}.npy'
    if f.exists():
        return np.load(f)
    t = tn.inference(probs).astype(np.float32)
    np.save(f, t)
    return t


def build_logits(seqstr, p16, tt, variant):
    """15-role fixed offset for one oriented strand (float32, shape [L,15])."""
    L = len(seqstr)
    p = np.clip(p16.astype(np.float64), MIN_PROB, MAX_PROB)
    d = (np.log(p) - np.log1p(-p)).astype(np.float32)
    row = np.zeros((L, 15), dtype=np.float32)
    for r in (1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14):
        row[:, r] = d
    if variant == 'B_emission':
        return row
    arr = np.frombuffer(seqstr.encode('ascii'), dtype=np.uint8)
    A, T, C, G = (arr == ord(c) for c in 'ATCG')
    # enter_gene: ATG starting at i
    atg = np.zeros(L, bool)
    atg[:-2] = A[:-2] & T[1:-1] & G[2:]
    row[atg, 13] += tt[atg, 0]
    # leave_gene at the last base j of an in-frame-capable stop triple
    stop = np.zeros(L, bool)
    for s in STOPS:
        m = np.ones(L - 2, bool)
        for k, ch in enumerate(s):
            m &= arr[k:L - 2 + k] == ord(ch)
        stop[2:] |= m
    row[stop, 14] += tt[stop, 1]
    # enter_intron: last exon base e carries the cost of the intron starting at e+1 (needs seq[e+1:e+3] = GT or GC)
    ent = np.zeros(L, np.float32)
    ent[:-1] = tt[1:, 2]
    ok_gt = np.zeros(L, bool); ok_gc = np.zeros(L, bool)
    ok_gt[:-2] = G[1:-1] & T[2:]
    ok_gc[:-2] = G[1:-1] & C[2:]
    ent[ok_gc] += np.float32(GC_AG_INTRON_LOG_PRIOR)
    ok = ok_gt | ok_gc
    for r in (7, 8, 9):
        row[ok, r] += ent[ok]
    # leave_intron: first exon base after the acceptor carries tr[i,3] (roles 10-12 occur only there)
    for r in (10, 11, 12):
        row[:, r] += tt[:, 3]
    return row


def run_strand(seq, strand, oriented, p16, tt, variant, max_exons=128):
    logits = torch.from_numpy(build_logits(oriented, p16, tt, variant))
    t0 = time.time()
    res = decode_full_grammar_skyline(oriented, logits, max_exons)
    return res, time.time() - t0


def carbon_region_keys(seq, sub_seq, tracks, tts, offset, cn_dev):
    """Carbon-A post-processing (numba Viterbi + CDS>=96) on exactly the same tracks/transitions; returns keys and genes."""
    import post_processing as PP
    from viterbi import choose_viterbi_backend
    from confidence_network import load_confidence_network
    codon = load_codon_table(str(CODON_TABLE_PATH))
    stops = tuple(stop_codons_from_table(codon))
    n = len(sub_seq)
    fw, rv = encode_both_strands(sub_seq)
    tr = [(tracks['+'][1], fw), (tracks['-'][1], rv)]
    ctx = [(seq, sub_seq, '+', n), (seq, sub_seq, '-', n)]
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, cn_dev)
    genes, short, vs, cs = PP.decode_chunk_genes(tr, ctx, [tts['+'], tts['-']], 'numba', stops, MAX_BOUNDARY_INTRON_LENGTH_BP,
                                                MIN_CDS_LENGTH_BP, cn)
    out = {}
    for g, _, _ in genes:
        k = (seq, g['strand'], tuple((s + offset, e + 1 + offset) for s, e in g['exons']))
        out[k] = dict(complete=g['complete'], reason=g.get('incomplete_reason'), conf=g['confidence'], exons=g['exons'])
    return out


def to_transcripts(seqid, sequence, strand, res, ordinal, offset=0):
    txs = []
    n = len(sequence)
    for i, (chain, gain) in enumerate(zip(res['chains'], res['chain_gains'])):
        g = genomic_chain(chain, 0, n, strand)
        g = [(a + offset, b + offset) for a, b in g]
        txs.append(as_transcript(g, seqid, strand, f'M29_GATE_{ordinal}_{"p" if strand == "+" else "m"}_{i}', gain))
    return txs


def main():
    variants = sys.argv[1].split(',')
    seqs = sys.argv[2].split(',')
    region = None
    if len(sys.argv) > 3 and sys.argv[3] != 'full':
        a, b = sys.argv[3].split('-')
        region = (int(a) * 1000, int(b) * 1000)    # in kb
    tn = load_transition_network(TRANSITION_NETWORK_PATH, torch.device('cpu'))
    for seq in seqs:
        rid, sequence, pos, neg = load_tracks(seq)
        assert rid == seq
        key = seq
        a, b = (0, len(sequence)) if region is None else region
        if region is not None:
            sub_seq, sub_pos, sub_neg = sequence[a:b], pos[a:b], neg[a:b]
            key = f'{seq}_{a}_{b}'
        else:
            sub_seq, sub_pos, sub_neg = sequence, pos, neg
        n = len(sub_seq)
        tracks = {'+': (sub_seq, sub_pos), '-': (rc(sub_seq), np.ascontiguousarray(sub_neg[::-1]))}
        tts = {s: transition_tracks(tn, tracks[s][1], f'{key}_{s}') for s in '+-'}
        for variant in variants:
            out = OUTROOT / variant / key
            out.mkdir(parents=True, exist_ok=True)
            txs, stats = [], {}
            for si, strand in enumerate('+-'):
                res, sec = run_strand(seq, strand, tracks[strand][0], tracks[strand][1], tts[strand], variant)
                txs += to_transcripts(seq, sub_seq, strand, res, 0, offset=a)
                stats[strand] = dict(chains=len(res['chains']), gain=res['gain'], seconds=sec)
                print(key, variant, strand, len(res['chains']), f'{sec:.0f}s', flush=True)
            gff = out / 'full_grammar.gff3'
            if gff.exists():
                gff.unlink()
            write_gff3(gff, txs, {seq: len(sequence)}, 'M29_GATE')
            (out / 'summary.json').write_text(json.dumps(dict(variant=variant, seq=seq, region=(a, b), chains=len(txs), stats=stats), indent=1))
            if region is not None:
                ca = carbon_region_keys(seq, sub_seq, tracks, tts, a, torch.device('cpu'))
                ours = {M.key(t) for t in txs}
                ca_complete = {k for k, v in ca.items() if v['complete']}
                only_ca = ca_complete - ours; only_ours = ours - ca_complete
                print(key, variant, 'CarbonA genes', len(ca), 'complete', len(ca_complete), 'ours', len(ours),
                      'identical', len(ours & ca_complete), 'only_CA', len(only_ca), 'only_ours', len(only_ours), flush=True)
                (out / 'parity_vs_carbon_same_tracks.json').write_text(json.dumps(dict(
                    carbon_genes=len(ca), carbon_complete=len(ca_complete), ours=len(ours), identical=len(ours & ca_complete),
                    only_carbon=[[list(map(list, k[2])), k[1]] for k in sorted(only_ca)][:200],
                    only_ours=[[list(map(list, k[2])), k[1]] for k in sorted(only_ours)][:200]), indent=0))


if __name__ == '__main__':
    main()
