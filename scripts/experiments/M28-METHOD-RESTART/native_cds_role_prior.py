"""C0 plus two frozen native per-base CDS margins, relative to background."""
import torch
from torch import nn
from torch.nn import functional as F
from src.m28.core import M28Core


@torch.no_grad()
def native_cds_log_odds(heads, token_features):
    if len(heads) != 2:
        raise ValueError('Both original native classification heads are required')
    batch, tokens, _ = token_features.shape
    scores = [head(token_features.to(torch.bfloat16)).reshape(batch, tokens*6, 2).float()
        for head in heads]
    return torch.stack([score[..., 1]-score[..., 0] for score in scores], dim=-1)


class NativeCDSRolePrior(M28Core):
    def __init__(self, feature_dim=2048, hidden=128):
        super().__init__('C0', feature_dim, hidden)
        # Background stays the unchanged reference. The other fourteen roles
        # learn separate coefficients for both native heads; no head selection.
        self.native_cds_role_prior = nn.Linear(2, 14, bias=False)
        nn.init.zeros_(self.native_cds_role_prior.weight)

    def load_c0_weights(self, state):
        missing = self.load_state_dict(state, strict=False)
        if missing.missing_keys != ['native_cds_role_prior.weight'] or missing.unexpected_keys:
            raise ValueError('Exactly the new zero prior may be absent from the original C0 state')

    def forward(self, token_features, dna, valid_bases, native_log_odds):
        if native_log_odds.shape != (*dna.shape[:2], 2):
            raise ValueError('Both native CDS margins must align to every original DNA base')
        output = super().forward(token_features, dna, valid_bases)
        prior = F.pad(self.native_cds_role_prior(native_log_odds), (1, 0))
        output['segmentation_logits'] = output['segmentation_logits'] + prior*output['mask'][..., None]
        return output
