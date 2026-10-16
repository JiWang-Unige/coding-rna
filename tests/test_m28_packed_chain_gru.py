"""CPU numerical equivalence only: no R4 data/checkpoints or throughput claim."""
import copy
import json
import pytest
import torch
from src.m28.core import ChainHead
from src.m28.packed_chain_gru import packed_chain_forward
from src.m28.training import balanced_binary_loss

torch.set_num_threads(2)
LONG = [(2+4*i, 4+4*i) for i in range(128)]
MIXED = [[(2, 11)], [(2, 11), (22, 34), (45, 61)],
         [(70, 79), (88, 99)], [(2, 11), (22, 34)]]
CASES = {
    "empty": ([], []),
    "single": ([[(3, 28)]], [1]),
    "mixed_unsorted_shared_exons": (MIXED, [0, 1, -1, 0]),
    "all_unknown": (MIXED, [-1, -1, -1, -1]),
    "long_128_exons": ([MIXED[0], LONG, MIXED[1]], [0, 1, -1]),
    "all_reference_negative": (MIXED, [0, 0, 0, 0]),
}


def check_case(name, chains, labels, dtype, hidden=16):
    # Fixed before the first run; never widen these after observing an error.
    atol, rtol = (2e-6, 2e-4) if dtype == torch.float32 else (1e-11, 1e-9)
    errors = {}
    def close(group, actual, expected):
        torch.testing.assert_close(actual, expected, atol=atol, rtol=rtol)
        if actual.numel():
            diff = (actual-expected).detach().abs()
            value = errors.setdefault(group, {"max_abs": 0., "max_relative_with_atol_floor": 0.})
            value["max_abs"] = max(value["max_abs"], float(diff.max()))
            value["max_relative_with_atol_floor"] = max(
                value["max_relative_with_atol_floor"],
                float((diff / expected.detach().abs().clamp_min(atol)).max()))
    torch.manual_seed(812)
    baseline = ChainHead(hidden).to(dtype)
    features = torch.randn(640, hidden, dtype=dtype)
    opt = torch.optim.AdamW(baseline.parameters(), lr=3e-4, weight_decay=.01,
                           betas=(.9, .999), eps=1e-8)
    # Populate nonzero Adam moments using only synthetic inputs and baseline code.
    warm = balanced_binary_loss(baseline(features, MIXED, 600)["logits"], [0, 1, -1, 0])
    warm.backward()
    torch.nn.utils.clip_grad_norm_(baseline.parameters(), 1.)
    opt.step()
    opt.zero_grad(set_to_none=True)
    candidate = copy.deepcopy(baseline)
    opt_new = torch.optim.AdamW(candidate.parameters(), lr=3e-4, weight_decay=.01,
                               betas=(.9, .999), eps=1e-8)
    opt_new.load_state_dict(copy.deepcopy(opt.state_dict()))
    before = {key: p.detach().clone() for key, p in baseline.named_parameters()}
    x = features.clone().requires_grad_(True)
    y = features.clone().requires_grad_(True)
    old = baseline(x, chains, 600)
    new = packed_chain_forward(candidate, y, chains, 600)
    for key in old:
        close("score", new[key], old[key])
    loss_old = 1.37 * balanced_binary_loss(old["logits"], labels)
    loss_new = 1.37 * balanced_binary_loss(new["logits"], labels)
    close("loss", loss_new, loss_old)
    loss_old.backward()
    loss_new.backward()
    close("feature_gradient", y.grad, x.grad)
    assert torch.count_nonzero(x.grad[600:]) == 0
    assert torch.count_nonzero(y.grad[600:]) == 0
    old_params, new_params = dict(baseline.named_parameters()), dict(candidate.named_parameters())
    for key, p in old_params.items():
        q = new_params[key]
        assert (p.grad is None) == (q.grad is None), key
        if p.grad is not None:
            close("parameter_gradient", q.grad, p.grad)
    norm_old = torch.nn.utils.clip_grad_norm_(baseline.parameters(), 1.)
    norm_new = torch.nn.utils.clip_grad_norm_(candidate.parameters(), 1.)
    close("gradient_norm", norm_new, norm_old)
    opt.step()
    opt_new.step()
    for key, p in old_params.items():
        q = new_params[key]
        close("parameter_increment", q.detach()-before[key], p.detach()-before[key])
        assert opt.state[p].keys() == opt_new.state[q].keys()
        for state_key in opt.state[p]:
            close("optimizer_state", opt_new.state[q][state_key], opt.state[p][state_key])
    print(json.dumps({"case": name, "dtype": str(dtype), "hidden": hidden,
                      "candidates": len(chains), "atol": atol, "rtol": rtol,
                      "positive_sign_flips": int(((old["logits"]>0)!=(new["logits"]>0)).sum()),
                      "errors": errors}))


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("dtype", [torch.float64, torch.float32])
def test_equivalence(name, dtype):
    chains, labels = CASES[name]
    check_case(name, chains, labels, dtype)


def test_actual_hidden_size_candidate_load():
    chains = [[(2+3*(j%40)+18*i, 11+3*(j%40)+18*i)
               for i in range(1+j%4)] for j in range(201)]
    check_case("201_candidates_hidden128", chains, [j%3-1 for j in range(201)], torch.float32, 128)


@pytest.mark.parametrize("dtype", [torch.float64, torch.float32])
def test_near_zero_gain_is_reported(dtype):
    torch.manual_seed(812)
    head = ChainHead(16).to(dtype)
    x = torch.randn(640, 16, dtype=dtype)
    with torch.no_grad():
        head.null.bias.add_(head(x, MIXED, 600)["logits"][0])
        old = head(x, MIXED, 600)["logits"]
        new = packed_chain_forward(head, x, MIXED, 600)["logits"]
    atol, rtol = (2e-6, 2e-4) if dtype == torch.float32 else (1e-11, 1e-9)
    torch.testing.assert_close(new, old, atol=atol, rtol=rtol)
    print(json.dumps({"case": "near_zero", "dtype": str(dtype),
                      "baseline_logits": old.tolist(), "packed_logits": new.tolist(),
                      "positive_sign_flips": int(((old>0)!=(new>0)).sum())}))
