"""Head/feature plumbing shared by training and inference: frozen NT-v2 features (selected hidden layers) -> M28 C0 head -> 15-role
potentials owned per base by the nearest NT-window centre (frozen `owned_context_logits`)."""
import re

import torch

import ntv2_adapter as A
import m31_common as M

G, Q = M.G, M.Q
# layers fed to the head; set by c1_train / c1_infer from the run config (hidden_states indices; 29 = last_hidden_state)
DEFAULT_LAYERS = (29,)


def feature_dim(layers):
    return A.HIDDEN * len(layers)


def new_head(layers, seed=0):
    return A.new_head(feature_dim(layers), seed)


def oriented_windows(chromosome, geometry, strand):
    """Oriented window strings of an NT geometry (genomic windows, increasing genomic start)."""
    out = []
    for left, right in geometry['genomic_windows']:
        s = re.sub('[^ACGT]', 'N', chromosome[left:right])
        out.append(Q.M.rc(s) if strand == '-' else s)
    return out


def window_potentials(model, feats, sequences, device='cuda'):
    """Valid-length 15-role potentials of several NT windows in ONE batched head call (windows are independent: the BiGRU runs on
    each window's own packed length), list of [len(s),15] tensors.  Autograd follows the caller's grad mode."""
    dna = torch.stack([G.one_hot(s, A.NT_WIDTH) for s in sequences]).to(device)
    f = torch.stack([x for x in feats]).float()
    out = model(f, dna, [len(s) for s in sequences])['segmentation_logits']
    return [out[i, :len(s)] for i, s in enumerate(sequences)]


def head_context_logits(model, feats, sequences, geometry, strand, device='cuda'):
    """Oriented-context potentials [len(context), 15] on `device`; autograd follows the grad mode of the caller."""
    return G.owned_context_logits(window_potentials(model, feats, sequences, device), geometry, strand)
