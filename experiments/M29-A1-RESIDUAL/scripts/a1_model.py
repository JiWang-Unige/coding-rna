"""Residual head: C0 (366,191 parameters, frozen-encoder features) whose role layer is zero-initialised, plus the fixed offset."""
import re
import torch
import a1_common as K
import fit_full_grammar_objective_contrast as C
import annotate_full_grammar_joint24 as A
K.repoint()
G, Q = C.G, C.Q


def new_model(seed=0):
    torch.manual_seed(seed)
    m = G.M28Core('C0', 2048, 128)
    torch.nn.init.zeros_(m.segmentation.weight)
    torch.nn.init.zeros_(m.segmentation.bias)
    assert sum(p.numel() for p in m.parameters()) == 366191
    return m


def oriented_windows(chromosome, geometry, strand):
    out = []
    for left, right in geometry['genomic_windows']:
        s = re.sub('[^ACGT]', 'N', chromosome[left:right])
        out.append(Q.M.rc(s) if strand == '-' else s)
    return out


@torch.inference_mode()
def encode_bf16(backbone, tokenizer, sequences, device='cuda'):
    """Frozen BF16 encoder on each 24,576-bp window; stored as BF16 (the values are already BF16-rounded, FP32 view is exact)."""
    return [A.encode_window(backbone, tokenizer, s, device).to(torch.bfloat16) for s in sequences]


def head_context_logits(model, feats, sequences, geometry, strand, device='cuda'):
    """Oriented-context head logits [len(context), 15] on `device` with autograd (if enabled), same ownership as the frozen path."""
    pieces = []
    for f, s in zip(feats, sequences):
        features = f.float().clone()
        dna = G.one_hot(s, Q.M.WIDTH)[None].to(device)
        pieces.append(model(features, dna, [len(s)])['segmentation_logits'][0, :len(s)])
    return G.owned_context_logits(pieces, geometry, strand)
