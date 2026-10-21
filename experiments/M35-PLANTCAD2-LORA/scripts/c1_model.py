"""Head/feature plumbing shared by training and inference: frozen PlantCAD2-Small per-base features (genomic orientation) -> oriented six-base
tokens -> M28 C0 head -> 15-role potentials owned per base by the nearest window centre (frozen `owned_context_logits`)."""
import re

import torch

import pc2_adapter as A
import m35_common as M

G, Q = M.G, M.Q


def new_head(seed=0):
    return A.new_head(A.HIDDEN, seed)


def genomic_window_strings(chromosome, geometry):
    """ACGTN strings of the genomic windows of a geometry (increasing genomic start)."""
    return [re.sub('[^ACGT]', 'N', chromosome[left:right]) for left, right in geometry['genomic_windows']]


def orient_strings(genomic, strand):
    return list(genomic) if strand == '+' else [Q.M.rc(s) for s in genomic]


def window_potentials(model, tokens, sequences, device='cuda'):
    """Valid-length 15-role potentials of several oriented windows in ONE batched head call (windows are independent: the BiGRU runs on
    each window's own packed length), list of [len(s),15] tensors.  Autograd follows the caller's grad mode."""
    dna = torch.stack([G.one_hot(s, A.HEAD_W) for s in sequences]).to(device)
    f = torch.stack([x for x in tokens]).float()
    out = model(f, dna, [len(s) for s in sequences])['segmentation_logits']
    return [out[i, :len(s)] for i, s in enumerate(sequences)]


def head_context_logits(model, tokens, sequences, geometry, strand, device='cuda'):
    """Oriented-context potentials [len(context), 15] on `device`; autograd follows the grad mode of the caller."""
    return G.owned_context_logits(window_potentials(model, tokens, sequences, device), geometry, strand)
