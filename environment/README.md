# environment

- `requirements-generanno.txt`：M28 参考管线的固定依赖（Baobab conda 环境 `generanno`）。
- `setup_codingrna_torch.sh`、`setup_annevo_env.sh`、`setup_helixer_container_m1.sbatch`、`job_watch.sh`：早期环境/监控脚本（从旧 `scripts/` 根目录移来，未改动）。

Baobab 上 sbatch 通过 `source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh; conda activate generanno`、`HF_HOME=/home/users/j/jwang/.cache/huggingface` 等硬编码路径启动，重组未改变这些。
