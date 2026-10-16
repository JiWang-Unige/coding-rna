"""New storage branch: exact arrays, Python arithmetic and full DP fixtures."""
import unittest

import numpy as np
import torch

import test_m28_emission_lattice_cached_metadata as C
from src.m28.emission_lattice import _Margins
from src.m28.emission_lattice_compact import _CompactMargins, decode_compact_emission_lattice


class CompactEmissionLatticeTests(C.CachedMetadataEmissionLatticeTests):
    def compare(self, sequence, potentials, cap, oracle=False):
        previous = super().compare(sequence, potentials, cap, oracle=oracle)
        actual = decode_compact_emission_lattice(sequence, potentials, cap)
        self.assertEqual(actual, previous)
        self.assertIs(type(actual['gain']), float)
        self.assertTrue(all(type(value) is float for value in actual['chain_gains']))
        return actual

    def test_arrays_and_scalar_operations_preserve_FP64_bits_and_storage_bound(self):
        generator = torch.Generator().manual_seed(79)
        random = torch.randn(257, 15, generator=generator)
        dynamic = random * (2. ** torch.tensor([-80., 0., 80.] * 86)[:257])[:, None]
        for logits in (torch.zeros(257, 15), random, dynamic):
            old, new = _Margins(logits), _CompactMargins(logits)
            for key in ('gain', 'prefix', 'coding'):
                expected = np.asarray(getattr(old, key), dtype=np.float64)
                actual = getattr(new, key)
                self.assertEqual(actual.dtype, np.dtype('float64'))
                self.assertTrue(np.array_equal(expected.view(np.uint64), actual.view(np.uint64)))
            self.assertEqual(sum(getattr(new, key).nbytes for key in ('gain', 'prefix', 'coding')),
                             264 * len(logits) + 144)
            for a in range(0, 248, 7):
                for size in (1, 2, 3, 4, 8):
                    b = a + size
                    for phase in range(3):
                        for first, terminal in ((False, False), (False, True),
                                                 (True, False), (True, True)):
                            if (first or terminal) and size < 3:
                                continue
                            expected = old.exon_gain(a, b, phase, first, terminal)
                            actual = new.exon_gain(a, b, phase, first, terminal)
                            self.assertIs(type(actual), float)
                            self.assertEqual(actual.hex(), expected.hex())


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
