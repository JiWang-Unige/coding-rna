"""Resumable-checkpoint helpers of the M34 trainer (no GPU code, unit-tested on CPU).  Same content as the M31 training-B checkpoint:
head, LoRA state, optimizer, RNG states, step, config, draw position.  `lora_state` is the dict of LoRA tensors, `copy_lora(state)` loads one in place."""
import numpy as np
import torch

import m34_train_common as T

M = T.M


def snapshot(head, lora_state, opt, step, cfg, draw_position):
    return dict(model=head.state_dict(), lora=lora_state, optimizer=opt.state_dict(), step=step, rng_torch=torch.get_rng_state(),
                rng_numpy=np.random.get_state(), config=dict(cfg), draw_position=draw_position)


def save_resumable(path, head, lora_state, opt, step, cfg, draw_position):
    st = snapshot(head, lora_state, opt, step, cfg, draw_position)
    M.atomic_save(st, path)
    return st


def save_model(path, st):
    """Model-only file (what inference loads): head + LoRA + config."""
    M.atomic_save(dict(model=st['model'], lora=st['lora'], step=st['step'], config=st['config']), path)


def restore(path, head, copy_lora, opt, cfg):
    """Load a resumable checkpoint into the live objects.  Returns the step.  Refuses a changed config."""
    st = torch.load(path, map_location='cpu', weights_only=False)
    assert st['config'] == dict(cfg), 'config changed since checkpoint'
    head.load_state_dict(st['model'])
    copy_lora(st['lora'])
    opt.load_state_dict(st['optimizer'])
    torch.set_rng_state(st['rng_torch'])
    np.random.set_state(st['rng_numpy'])
    return st['step']
