"""M40-SEMISUP-MLM shared definitions (protocol docs/review_20261016/32_semisup_MLM_protocol.md).

Semi-supervised variant of the M34 pc2 trainer: labeled positions keep the partial-observation CRF likelihood (unknown positions marginalised, treatment U),
and the UNLABELED positions of the NEW species' draws additionally carry a masked-language-model loss through PlantCAD2's own pre-trained (tied, frozen)
RCPS LM head.  At label fraction 0 (arm M40_k0, treatment X) the new species carry no CRF term at all, only the MLM.

This module only IMPORTS the M34 code (M34 files are not edited): it sets the backend to pc2, puts the M34 scripts on the path and registers the M40 arms in
the in-process `m34_schedule.ARMS` dict, so that the M34 helpers (labelers, validation regions, `arm_dir`, `m34_select.py`) work for them unchanged.
"""
import hashlib
import os
import sys
from pathlib import Path

os.environ.setdefault('M34_BACKEND', 'pc2')
HERE = Path(__file__).resolve().parent
EXPS = HERE.parents[1]
M34_SCRIPTS = EXPS / 'M34-NTV2-SCALEUP' / 'scripts'
for _p in (str(M34_SCRIPTS), str(HERE)):
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)                         # HERE first, then the M34 scripts

import m34_train_common as T                       # noqa: E402  (pc2: M32 adapter / c1_model, M34 labels)
import m34_schedule as S                           # noqa: E402
import m34_common as C                             # noqa: E402

assert T.BACKEND == 'pc2', 'M40 is defined for the PlantCAD2 backend only (M34_BACKEND=pc2)'

OUT40 = C.OUT.parent / 'M40-SEMISUP-MLM'      # ~/coding-rna/outputs/M40-SEMISUP-MLM (repeat index, k0 schedule, smoke)
REPEATS = OUT40 / 'repeats'
SCHEDULES40 = OUT40 / 'schedules'

# MLM configuration (protocol 32 section 3; fixed before any M40 training output).  Vipsania reference values in brackets.
MLM_CFG = dict(
    mlm_rate=0.05,              # Bernoulli masking rate per ACGT base of a new-species draw context  [Vipsania: at most 5 % of the bases]
    mlm_mask_frac=0.8,          # of the selected bases: 80 % -> [MASK], 10 % -> random base, 10 % unchanged  [PlantCAD2 pre-training 80/10/10]
    mlm_rand_frac=0.1,
    mlm_repeat_weight=1e-4,     # loss weight of soft-masked (repeat) bases  [Vipsania: repeats down-weighted to 1e-4]
    mlm_lambda_max=0.1,         # weight of the MLM term relative to the per-known-base CRF NLL of a draw
    mlm_start_frac=0.1,         # lambda = 0 (and no masking) for the first 10 % of the steps  [Vipsania: auxiliary weight 0 first, then raised]
    mlm_ramp_frac=0.1,          # then linear ramp to lambda_max over the next 10 %, constant afterwards
    mlm_scope='unlabeled',      # loss only at positions without a CRF label (unknown core positions and the context flanks); masking itself is label-blind
    mlm_species='new',          # MLM only on draws of the 8 new species (the two legacy species are fully labeled and keep the M34 path)
    mlm_seed=40,                # per-draw mask RNG = default_rng([mlm_seed, draw index]) (resume-invariant)
)

LEGACY, NEW = list(S.LEGACY), list(S.NEW)
ARMS40 = {
    # k0: no label on any new-species base (treatment X); legacy species fully labeled as in every M34 arm.  Own draw schedule: the M34 sampler at k = 0
    # (the gene-first half of the sampler uses labeled genes, so a k = 0 arm cannot reuse a labeled arm's schedule without leaking labels into sampling).
    'M40_k0': dict(species=LEGACY + NEW, new_k=0.0, treatment='X', mlm=True, schedule='M40', control=None),
    # partial-label arms: the M34 arm (same labeled genes, same draw schedule file, same seeds / init) + the MLM term
    'M40_k25U_MLM': dict(species=LEGACY + NEW, new_k=0.25, treatment='U', mlm=True, schedule='T10_k25U', control='T10_k25U'),
    'M40_k50U_MLM': dict(species=LEGACY + NEW, new_k=0.5, treatment='U', mlm=True, schedule='T10_k50U', control='T10_k50U'),
}


def register():
    """Register the M40 arms in m34_schedule.ARMS (in-process only; the M34 labelers / val-region helpers / arm_dir then accept them)."""
    for k, v in ARMS40.items():
        S.ARMS.setdefault(k, {kk: vv for kk, vv in v.items() if kk in ('species', 'new_k', 'treatment')})
    return S.ARMS


register()


def schedule_path(arm):
    src = ARMS40[arm]['schedule']
    return SCHEDULES40 / f'schedule_{arm}.jsonl' if src == 'M40' else T.SCHEDULES / f'schedule_{src}.jsonl'


def md5_file(path):
    return hashlib.md5(open(path, 'rb').read()).hexdigest()


def species_treatment(arm, species):
    """Treatment used for the CRF targets of one species in an arm: legacy species are always fully labeled ('U' at k = 1 is identical to FULL)."""
    t = ARMS40[arm]['treatment']
    if species in LEGACY:
        return 'U' if t == 'X' else t
    return t


def mlm_on_species(species, cfg=MLM_CFG):
    return cfg['mlm_species'] == 'all' or species not in LEGACY
