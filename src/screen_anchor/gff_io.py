"""Convert per-base predicted labels -> a CDS GFF that scripts/eval_gene_body_mask.py
can read under --span-mode cds. Pure stdlib + numpy.

Grouping: a predicted "gene" = a maximal run of label>0 (CDS or gene-body-nc). Within it,
CDS sub-runs (label==1) are emitted as CDS lines sharing one transcript_id/gene_id, so the
evaluator's CDS-span for that gene = first CDS start .. last CDS end. Emitted on '+' strand
(models are strand-agnostic; the cds-span evaluator groups by (seqid, strand, group_id)).
"""
import numpy as np

from .data import CLASS_CDS


def _runs(mask):
    """Yield (start, end) 0-based half-open runs where boolean mask is True."""
    if mask.size == 0:
        return
    idx = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    for i in range(0, len(idx), 2):
        yield int(idx[i]), int(idx[i + 1])


def labels_to_cds_gff(pred_by_seqid, out_path, source="screen_ref"):
    """pred_by_seqid: {seqid: int8 array of predicted classes}. Writes GTF-attr CDS lines."""
    n_genes = 0
    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for seqid in sorted(pred_by_seqid):
            arr = pred_by_seqid[seqid]
            genebody = arr > 0
            cds = arr == CLASS_CDS
            for g0, g1 in _runs(genebody):
                cds_sub = cds[g0:g1]
                if not cds_sub.any():
                    continue  # no CDS in this gene-body run -> nothing to score in cds mode
                n_genes += 1
                gid = f"{seqid}_g{n_genes}"
                for c0, c1 in _runs(cds_sub):
                    start = g0 + c0 + 1            # GFF 1-based inclusive
                    end = g0 + c1
                    fh.write(
                        f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t+\t0\t"
                        f'transcript_id "{gid}"; gene_id "{gid}";\n'
                    )
    return n_genes
