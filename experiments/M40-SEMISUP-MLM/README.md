# M40-SEMISUP-MLM

Partial-label CRF (M34 pc2 trainer, unknown positions marginalised) + masked-language-model auxiliary loss on unlabeled new-species positions, through
PlantCAD2-Small's own pre-trained tied RCPS LM head.  Protocol and pre-registered criteria: `docs/review_20261016/32_semisup_MLM_protocol.md`.

- `scripts/m40_common.py` arms (M40_k0, M40_k25U_MLM, M40_k50U_MLM) and MLM config; `m40_mlm.py` masking / LM head / lambda schedule;
  `m40_repeats.py` soft-mask index; `m40_schedule.py` k0 schedule; `m40_train.py` trainer.
- Shell (ont): `m40_smoke.sh` (GPU smoke), `m40_stage.sh <arm>` (= m34_stage with the M40 trainer), `m40_run_all.sh` (queue entry: smoke -> arms).
- CPU prep (done 2026-11-11): `m40_schedule.py`; `m40_repeats.py build`.  Outputs: `~/coding-rna/outputs/M40-SEMISUP-MLM/`; arms in
  `~/coding-rna/outputs/M34-NTV2-SCALEUP/pc2_<arm>/`.
- Tests (ont CPU): `M34_BACKEND=pc2 PYTHONPATH=~/envs/m32_extra ~/envs/plantcad2/bin/python -m pytest -q tests/`.
- Push code: `./push.sh`.
