"""Scoped precision candidate motivated by D1; never imported by R4 runners."""
import torch
from .packed_chain_gru import packed_chain_forward


def packed_chain_forward_no_tf32(head, features, chains, valid_bases=None):
    # This single-threaded research process changes only the chain-call flag.
    # Do not use cudnn.flags(allow_tf32=False): its other defaults disable cuDNN.
    previous = torch.backends.cudnn.allow_tf32
    torch.backends.cudnn.allow_tf32 = False
    try:
        return packed_chain_forward(head, features, chains, valid_bases)
    finally:
        torch.backends.cudnn.allow_tf32 = previous
