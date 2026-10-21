import sys
import m31_common as M, a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_common as MC
for seq, a, b in (('NC_003071.7', 5_000_000, 13_000_000), ('NC_089044.1', 14_000_000, 26_000_000), ('NC_003071.7', 0, 10**9), ('NC_089044.1', 0, 10**9)):
    ref = MC.build_reference(seq)
    n = sum(1 for t in ref['refs'] if all(x >= a and y <= b for x, y, _ in t['CDS']))
    print(seq, a, b, 'refs_in_region', n, 'of', len(ref['refs']), 'length', MC.LENGTH[seq])
