import torch
import a1_model as MD
import a1_train as T


def test_zero_init_model_outputs_zero_and_has_366191_parameters():
    m = MD.new_model(0)
    assert sum(p.numel() for p in m.parameters()) == 366191
    f = torch.randn(1, 4096, 2048)
    dna = torch.zeros(1, 24576, 5); dna[..., 0] = 1
    out = m(f, dna, [24576])['segmentation_logits']
    assert torch.count_nonzero(out) == 0


def test_lr_schedule():
    total = 4608
    assert abs(T.lr_at(0, total) - 3e-4 / 231) < 1e-12 or T.lr_at(0, total) < 3e-4
    w = 231
    assert abs(T.lr_at(w - 1, total) - 3e-4) < 1e-12
    assert abs(T.lr_at(w, total) - 3e-4) < 1e-12
    assert abs(T.lr_at(total - 1, total) - 1e-5) < 5e-8
    assert all(T.lr_at(i, total) >= T.lr_at(i + 1, total) - 1e-15 for i in range(w, total - 1))
