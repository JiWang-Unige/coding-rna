"""Build the frozen-ruler reference objects (M26 CDS-assessable primaries) for the B contig and for A under the M37 cache; print counts."""
import json, sys
from collections import Counter
import m37_common as C
M = C.register()
import semantic_reference_pair as R7
out = {}
for seq in (C.B_CONTIG, C.A_CONTIG):
    ref = M.build_reference(seq)
    nex = Counter(C.strat(len(t['CDS'])) for t in ref['refs'])
    ann = ref['annotation']
    out[seq] = dict(primaries=len(ref['refs']), assessable_isoforms=len(ref['isoforms']), gene_spans=len(ref['gene_spans']),
                    annotated_genes=len(ann['genes']), annotated_transcripts=len(ann['transcripts']),
                    phase_consistent=all(R7.phase_chain_consistent(t) for t in ref['refs']), strata=dict(nex))
    print(seq, json.dumps(out[seq]), flush=True)
C.OUT.mkdir(parents=True, exist_ok=True)
json.dump(out, open(C.OUT / 'reference_counts.json', 'w'), indent=1)
