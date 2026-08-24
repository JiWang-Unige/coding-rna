"""M1-SAMEBUDGET-SCREEN-ANCHOR — two light random-init per-base classifiers.

Both map one-hot DNA (B, C=5, W) -> per-base logits (B, W, NUM_CLASSES=3).
Mechanism delta (orthogonality):
  - TiberiusLike : Conv stack -> biLSTM  (recurrent context; Tiberius-family inductive bias)
  - HelixerLike  : stacked DILATED convs, NO recurrence (wide receptive field via dilation;
                   Helixer-family segmentation inductive bias)
Kept small (same-budget screen): a few layers, modest width.
"""
import torch
import torch.nn as nn

from .data import NUM_CHANNELS, NUM_CLASSES


class TiberiusLike(nn.Module):
    def __init__(self, hidden=128, conv_channels=64, lstm_layers=2, dropout=0.1):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(NUM_CHANNELS, conv_channels, kernel_size=9, padding=4),
            nn.ReLU(),
            nn.Conv1d(conv_channels, conv_channels, kernel_size=9, padding=4),
            nn.ReLU(),
        )
        self.lstm = nn.LSTM(conv_channels, hidden, num_layers=lstm_layers,
                            batch_first=True, bidirectional=True, dropout=dropout)
        self.head = nn.Linear(2 * hidden, NUM_CLASSES)

    def forward(self, x):                 # x: (B, C, W)
        h = self.conv(x)                  # (B, conv_channels, W)
        h = h.transpose(1, 2)             # (B, W, conv_channels)
        h, _ = self.lstm(h)               # (B, W, 2*hidden)
        return self.head(h)               # (B, W, NUM_CLASSES)


class HelixerLike(nn.Module):
    def __init__(self, channels=96, n_blocks=6, dropout=0.1):
        super().__init__()
        layers = [nn.Conv1d(NUM_CHANNELS, channels, kernel_size=7, padding=3), nn.ReLU()]
        dilation = 1
        for _ in range(n_blocks):
            layers += [
                nn.Conv1d(channels, channels, kernel_size=3, padding=dilation, dilation=dilation),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            dilation = min(dilation * 2, 64)   # exponentially growing receptive field
        self.body = nn.Sequential(*layers)
        self.head = nn.Conv1d(channels, NUM_CLASSES, kernel_size=1)

    def forward(self, x):                 # x: (B, C, W)
        h = self.body(x)                  # (B, channels, W)
        return self.head(h).transpose(1, 2)  # (B, W, NUM_CLASSES)


def build_model(name, **kw):
    name = name.lower()
    if name in ("tiberius_like", "tiberius-like", "tiberiuslike"):
        return TiberiusLike(**kw)
    if name in ("helixer_like", "helixer-like", "helixerlike"):
        return HelixerLike(**kw)
    raise ValueError(f"unknown model '{name}'")
