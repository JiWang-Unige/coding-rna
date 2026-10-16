# tests：仅针对 `pipeline/m28` 内代码的测试

保留在这里的测试，其被测模块（`src.m28.*`、`src.screen_anchor.*`、`src.foundation_probe.*` 或 `pipeline/m28/scripts/*.py`）都位于 [`pipeline/m28`](../pipeline/m28/README.md)，
并且所有被它们 import/加载的文件也在 `pipeline/m28` 或本目录内（含 `test_m28_emission_lattice.py` 这类被其它测试 import 的共享夹具）。

其余测试随所测代码归档：

| 测试 | 位置 |
|---|---|
| M25R / M26（非流水线模块）/ M27 | `experiments/<group>/tests/` |
| M28 full-grammar 分支（LoRA、dense-QV、E/M、train-recovery、qv-lr 等） | [`experiments/M28/full_grammar/tests/`](../experiments/M28/full_grammar/tests/) |
| 基线运行器（native） | [`benchmarks/native_baselines/tests/`](../benchmarks/native_baselines/tests/) |
| A5–A44 及 pre-full-grammar | [`archive/m28_pre_full_grammar/tests/`](../archive/m28_pre_full_grammar/tests/) |

## 运行

```bash
# 仓库根；需要 torch/numpy（本机 Python 3.9 没有 torch，实际在 Baobab 的 generanno 环境运行）
export PYTHONPATH="$PWD/pipeline/m28:$PWD/pipeline/m28/scripts:$PWD/tests"
python -m pytest -q tests                                   # pyproject.toml 中已配置 pythonpath
python -m unittest discover -s tests -p 'test_m28_full_grammar_skyline.py' -v   # 与 sbatch 中的 checks 作业一致
```

注意：部分测试 `from test_m28_emission_lattice import ...` 之类跨测试 import，需要 `tests/` 在 `sys.path` 上（上面的 PYTHONPATH 或 pytest 配置提供）。
