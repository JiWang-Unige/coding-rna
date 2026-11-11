"""M41 inference = the M34 inference (`m34_infer.py`; with --setb `m34_p2b_infer.py`, which makes the evaluation-set-B contigs readable), imported unchanged,
with ONE change: every decode worker uses the M41 skyline decoder (`m41_grammar.decode_skyline`, short / long intron states, edge bonus of the model
file) instead of the frozen `m29_skyline_reuse.decode_full_grammar_skyline`.  Encoder, stitching, potentials, GFF3 export and file layout are the M34 ones,
so the M34 scorers (m34_val.py, c1_eval.py / c1_boot.py, m34_p2b_score.py) read the output as they read T2's.

The model file's `intron` entry (written by m41_train.py) is passed to the spawned decode workers through the environment variable M41_INTRON.
A model file without `intron` is refused (use m34_infer.py for frozen-grammar models).
Usage: M34_BACKEND=pc2 m41_infer.py [--setb] --model model_N.pt --regions A,B,... --out DIR [--decode-workers 4] [--cpus ...]
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T              # noqa: E402,F401  (backend first, as m34_p2b_infer requires)
SETB = '--setb' in sys.argv or os.environ.get('M41_SETB') == '1'      # spawned decode workers re-import this file: they see the env, not the flag
if SETB:
    import m34_p2b_infer as _B            # noqa: E402  (patches T.register for the M30 contigs, imports m34_infer)
    X = _B.X
else:
    import m34_infer as X                 # noqa: E402
if str(HERE) not in sys.path:
    sys.path.insert(1, str(HERE))

_ORIG_INIT, _ORIG_JOB = X._decode_init, X._decode_job


def _install():
    """Point the frozen decoder name used inside the M31/M32 `c1_infer._decode_job` at the M41 decoder (this process only)."""
    import torch
    import m29_skyline_reuse as SK
    import m41_grammar as GR
    if getattr(SK, '_m41_installed', False):
        return
    cfg = json.loads(os.environ['M41_INTRON'])
    il = GR.IntronLength(True, int(cfg['short_stages']))
    bonus = torch.tensor(cfg['bonus'], dtype=torch.float64)

    def decode(sequence, logits, max_exons, checkpoint_bp=8192):
        return GR.decode_skyline(sequence, logits, max_exons, bonus, il, checkpoint_bp)
    SK.decode_full_grammar_skyline = decode
    SK._m41_installed = True
    print(json.dumps(dict(event='m41_decoder_installed', pid=os.getpid(), short_stages=il.short_stages, state_count=len(GR.state_table(128, il.key())[0]))),
          file=sys.stderr, flush=True)


def _decode_init(cpus):
    _ORIG_INIT(cpus)
    _install()


def _decode_job(*args):
    _install()
    return _ORIG_JOB(*args)


X._decode_init, X._decode_job = _decode_init, _decode_job


def main():
    import torch
    if '--setb' in sys.argv:
        sys.argv.remove('--setb')
        os.environ['M41_SETB'] = '1'
    model = sys.argv[sys.argv.index('--model') + 1]
    st = torch.load(model, map_location='cpu', weights_only=False)
    if 'intron' not in st:
        raise SystemExit(f'{model} has no M41 intron-length entry; use m34_infer.py for frozen-grammar models')
    os.environ['M41_INTRON'] = json.dumps(dict(short_stages=st['intron']['short_stages'], bonus=st['intron']['bonus']))
    del st
    X.main()


if __name__ == '__main__':
    main()
