"""R1: Carbon-A transition-network ablation on saved raw probabilities (CPU, no re-inference).

Uses an UNMODIFIED COPY of annotation_pipeline/src in  experiments/M29-CARBON-A-GAP/carbon_a_pipeline_copy/ .
Raw probabilities: outputs/M28-CARBON-A-TRAIN-RAW/<seq>/native/*.raw_probs.parquet  (the saved GPU/bf16 inference).
For each variant the four per-base transition log-odds tracks are replaced (the transition TCN is bypassed) and the
rest of post-processing (30-state Viterbi, confidence TCN, CDS>=96, record>=100 kb) runs unchanged on CPU/fp32.

Variants
  control          : transition TCN as shipped (CPU fp32 replay; differs from saved GPU/bf16 GFF only by numerics)
  g<G>_i<I>        : position-independent costs. G = total cost per gene (half on enter_gene, half on leave_gene),
                     I = total cost per intron (half on enter_intron, half on leave_intron); either may be 'net' = keep the TCN track.
                     g0_i0 is the all-zero ablation. Costs are in nats of log-odds (same units as the TCN output, which spans about +-200).
Output: outputs/M29-CARBON-A-GAP/R1/<variant>/<seq>/prediction.gff  and  summary.json
"""
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = HERE.parent / 'carbon_a_pipeline_copy' / 'src'
sys.path.insert(0, str(SRC))
import torch  # noqa: E402

NT = int(os.environ.get('M29_THREADS', '5'))
torch.set_num_threads(NT)
os.environ.setdefault('OMP_NUM_THREADS', str(NT))
from constants import (CHUNK_SIZE_MBP, CODON_TABLE_PATH, CONFIDENCE_NETWORK_PATH, MAX_BOUNDARY_INTRON_LENGTH_BP,  # noqa: E402
                       MIN_CDS_LENGTH_BP, MIN_RECORD_LENGTH_BP, NEGATIVE_STRAND_PROBABILITY_COLUMN,
                       POSITIVE_STRAND_PROBABILITY_COLUMN, RECORD_ID_COLUMN, SEQUENCE_COLUMN, TRANSITION_NETWORK_PATH)
from utils import load_codon_table, stop_codons_from_table, write_gff3  # noqa: E402
from confidence_network import load_confidence_network  # noqa: E402
from transition_network import load_transition_network  # noqa: E402
from viterbi import choose_viterbi_backend  # noqa: E402
import post_processing as PP  # noqa: E402

GRID_G = (0, -20, -60, -150)
GRID_I = (0, -10, -30, -80)
VARIANTS = ['control'] + [f'g{g}_i{i}' for g in GRID_G for i in GRID_I] + [f'g{g}_inet' for g in GRID_G] + [f'gnet_i{i}' for i in GRID_I]
OUTROOT = ROOT / 'outputs/M29-CARBON-A-GAP/R1'
RAW = ROOT / 'outputs/M28-CARBON-A-TRAIN-RAW'
EV = {0: 'enter_gene', 1: 'leave_gene', 2: 'enter_intron', 3: 'leave_intron'}


class DummyConfidence:
    """Skips the confidence TCN (slow); every complete gene gets confidence 1.0, so 'filtered' == 'unfiltered' for these variants."""
    def inference(self, probabilities, exons, scored):
        return [1.0] * len(scored)


def realized_means(chunk_genes, tt_by_prob_id):
    """mean track value at the events used by the control decode (oriented coordinates)."""
    acc = {k: [] for k in range(4)}
    for gene, probs, exons in chunk_genes:
        tt = tt_by_prob_id[id(probs)]
        acc[0].append(tt[exons[0][0], 0])
        acc[1].append(tt[exons[-1][1], 1])
        for (a, b), (c, d) in zip(exons, exons[1:]):
            acc[2].append(tt[b + 1, 2])
            acc[3].append(tt[c, 3])
    return {k: float(np.mean(v)) if v else 0.0 for k, v in acc.items()}, {k: len(v) for k, v in acc.items()}


def apply_variant(variant, tt, consts=None):
    if variant == 'control':
        return tt
    g, i = variant[1:].split('_', 1)[0], variant.split('_', 1)[1][1:]
    out = tt.copy()
    if g != 'net':
        out[:, 0] = float(g) / 2; out[:, 1] = float(g) / 2
    if i != 'net':
        out[:, 2] = float(i) / 2; out[:, 3] = float(i) / 2
    return out


def run_contig(seq, variants, backend, tn, cn, codon):
    t0 = time.time()
    pq_path = sorted((RAW / seq / 'native').glob('*.raw_probs.parquet'))[0]
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(pq_path)
    stops = tuple(stop_codons_from_table(codon))
    chunks = list(PP.iter_probability_chunks(pf, RECORD_ID_COLUMN, SEQUENCE_COLUMN, POSITIVE_STRAND_PROBABILITY_COLUMN,
                                              NEGATIVE_STRAND_PROBABILITY_COLUMN, CHUNK_SIZE_MBP * 1_000_000,
                                              min_record_length=MIN_RECORD_LENGTH_BP, filter_stats={'record_count': 0, 'bp': 0}))
    assert len(chunks) == 1, len(chunks)
    tracks, contexts, bp = chunks[0]
    cache = ROOT / 'outputs/M29-CARBON-A-GAP/_cache'
    cache.mkdir(parents=True, exist_ok=True)
    tt = []
    for strand, (p, _) in zip('+-', tracks):
        f = cache / f'tt_{seq}_{strand}.npy'          # shared with m29_gate_parity.py
        if f.exists():
            tt.append(np.load(f))
        else:
            x = tn.inference(p).astype(np.float32)
            np.save(f, x)
            tt.append(x)
    print(seq, 'transition tracks', round(time.time() - t0, 1), 's', flush=True)
    consts, counts = None, None
    summary = {}
    for v in variants:
        t1 = time.time()
        out = OUTROOT / v / seq
        if v != 'control' and (out / 'prediction.gff').exists() and not os.environ.get('M29_FORCE'):
            print(seq, v, 'exists, skipped', flush=True)
            continue
        if v == 'control' and (out / 'prediction.gff').exists():
            print(seq, v, 'exists, skipped (control)', flush=True)
            continue
        use = [apply_variant(v, x) for x in tt] if v != 'control' else tt
        cn_use = cn if (v == 'control' or os.environ.get('M29_CONF')) else DummyConfidence()
        genes_all, short, vs, cs = PP.decode_chunk_genes(tracks, contexts, use, backend, stops, MAX_BOUNDARY_INTRON_LENGTH_BP,
                                                          MIN_CDS_LENGTH_BP, cn_use)
        if v == 'control' and consts is None:
            consts, counts = realized_means(genes_all, {id(p): x for (p, _), x in zip(tracks, tt)})
            summary['control_track_value_at_decoded_events_mean'] = consts; summary['decoded_event_counts'] = counts
            summary['track_mean_all_bases'] = {EV[k]: float(np.mean([x[:, k].mean() for x in tt])) for k in range(4)}
            summary['track_std_all_bases'] = {EV[k]: float(np.mean([x[:, k].std() for x in tt])) for k in range(4)}
        genes = sorted((g for g, _, _ in genes_all), key=lambda g: (g['record_id'], g['exons'][0][0]))
        out = OUTROOT / v / seq
        out.mkdir(parents=True, exist_ok=True)
        write_gff3(genes, str(out / 'prediction.gff'))
        summary[v] = dict(genes=len(genes), short_omitted=short, viterbi_s=vs, conf_s=cs, wall_s=time.time() - t1)
        print(seq, v, 'genes', len(genes), f'viterbi {vs:.0f}s conf {cs:.0f}s', flush=True)
    (OUTROOT / f'_run_{seq}_{os.getpid()}.json').write_text(json.dumps(summary, indent=1))


def main():
    variants = sys.argv[1].split(',') if sys.argv[1] != 'all' else VARIANTS
    backend = choose_viterbi_backend('cpu')
    print('backend', backend, 'threads', NT, flush=True)
    codon = load_codon_table(str(CODON_TABLE_PATH))
    dev = torch.device('cpu')
    tn = load_transition_network(TRANSITION_NETWORK_PATH, dev)
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, dev)
    for seq in sys.argv[2:]:
        run_contig(seq, variants, backend, tn, cn, codon)


if __name__ == '__main__':
    main()
