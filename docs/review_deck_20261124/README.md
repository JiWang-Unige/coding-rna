# 系统梳理幻灯片 第二版（2026-11-24）

这个文件夹放 2026-11-24 的系统梳理幻灯片及配套文档，供研究者本人捋思路、和 Claude 一起定论文故事用。它是 `docs/review_deck_20261102/`（11-02，28 页）的增量更新。
- 11-02 版文件夹保持不变，作为历史快照。
- 研究日志和各实验的协议、报告仍写在 `docs/review_20261016/`。

- **幻灯片**：https://claude.ai/artifact/Vfna3uyDGZz19dSz6s4N7m
  - 与 11-02 版是同一个链接，现为 36 页，版本 15。
  - Slides artifact，私有，别人要等你从页面的 Share 菜单分享后才能打开。

## 和 11-02 版相比

- **新增 12 页**：
  - 项目历程（journey）与遇到的问题（problems）；
  - C1 面板判定（c1-panel）；
  - P4 两 seed（c2-p4，重画）与训练动态（c2-dyn）；
  - 外部审阅的两点批评（c2-critique）；
  - M42 设计与关卡 0（c2-m42design）、关卡 1（c2-m42）、合并分析（c2-merge）；
  - C3 终版的可靠性（c3-rel）与选择性注释（c3-sel）；
  - 长基因（c4-long）与 M41（c4-m41）。
- **重写**：封面、变化、主张图、已放弃方向、综合、状态、决策分支、时间线、局限、待讨论问题。
- **删除**：c1-scale（并入"已放弃方向"）、c2-ckpt（由 c2-dyn 取代）、c3-pc2、branches2。
- **原 c1-panel 改名为 c1-panelset**：面板物种构成。
- **其余保留页**：只修正过时句子。修正内容见 `scripts/build_1124.py` 里的 `patch()` 调用。

## 内容

| 路径 | 内容 |
|---|---|
| `coding-rna_系统梳理_20261124.pdf` | 36 页幻灯片 PDF |
| `slides_png/` | 每页一张 PNG |
| `deck_source/` | 幻灯片源文件（`project/deck.json`、`project/slides/*.html`，含英文讲稿）和 `offline.html` |
| `幻灯片逐页详解.md` / `.pdf` | 中文逐页详解，只写当前状态 |
| `讲解稿_系统梳理.md` | 中文讲解稿，顶部是带日期的更新日志，保留 11-02 版的全部日志 |
| `verification/` | 四个只读核实的明细（A C1 面板与 M41、B C2、C C3 与 C4，TSV 每行带 path:line）和项目历程整理（D_history.md） |
| `figures/` | 18 张图：11-02 版的 11 张重新生成，另加 7 张新图 |
| `scripts/` | `style.py`（含 `cite()` / `has()` / `find()` 出处断言）、每张图一个 `fig_*.py`、`build_1124.py`（由 11-02 版生成本版的全部改动） |
| `blobs.json`、`blobs_flat.json` | 已上传图片的 asset id 到本地文件的对应（含 sha256） |

## 重建

1. 图（需要 matplotlib，不连任何远程机器）：

   ```bash
   python docs/review_deck_20261124/scripts/fig_p4dyn.py
   ```

2. 幻灯片：把 11-02 版的 `deck_source/project` 复制到一个临时目录，运行 `python scripts/build_1124.py <该目录>`，再运行 skill 的 `renumber.py` 和 `lint_slides.py`。
3. 导出：

   ```bash
   python ~/.claude/skills/research-review-deck/scripts/export_deck.py --deck <project> --out docs/review_deck_20261124 --name "coding-rna_系统梳理_20261124" --blobs docs/review_deck_20261124/blobs.json
   ```
