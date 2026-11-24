"""fig_p4seeds  (SLOT_WIDE 1664x560)  slide c2-p4

Claim C2: with half of the new species' genes labelled, marginalising (k50U) beats fill-background
(k50F) on all 8 new species in both seeds; seed 2's margin is mostly k50F_s2 under-calling.
Source: docs/review_20261016/20_M34_stage1_results.md §16.2 table (per-species s2 diff, s1 diff, mean).
"""
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R20 = "docs/review_20261016/20_M34_stage1_results.md"
first = style.find(R20, "| Populus trichocarpa | 72.24")
rows = []
for ln in range(first, first + 9):
    c = style.cells(R20, ln)
    u2, f2, d2, d1, m = (style._num(c[i]) for i in (1, 2, 3, 4, 5))
    style.close(u2 - f2, d2, 0.011)
    style.close((d1 + d2) / 2, m, 0.011)
    rows.append((c[0], d1, d2, m))
pooled = rows.pop()
assert pooled[0] == "pooled" and len(rows) == 8
style.has(R20, style.find(R20, "**判定：成立**"), "+9.10")
rows.sort(key=lambda r: -r[3])

W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
ax = style.axes_px(fig, 330, 70, 1150, 400)
data = rows + [pooled]
ys = list(range(len(rows) + 1, 1, -1)) + [0]
for y, (name, d1, d2, m) in zip(ys, data):
    bold = name == "pooled"
    ax.barh(y + 0.18, d1, height=0.34, color=style.BLUE_LIGHT)
    ax.barh(y - 0.18, d2, height=0.34, color=style.BLUE)
    ax.scatter([m], [y], marker="D", s=70, color=style.INK, zorder=3)
    ax.text(-0.8, y, "Pooled (8 species)" if bold else name, ha="right", va="center", fontsize=18,
            fontstyle="normal" if bold else "italic", fontweight="bold" if bold else "normal")
    ax.text(27.5, y, f"{d1:+.2f} / {d2:+.2f}", ha="left", va="center", fontsize=17,
            fontweight="bold" if bold else "normal", color=style.BODY)
ax.axvline(0, color=style.MUTED, lw=1.2)
ax.axhline(1.0, color=style.GRID, lw=1.2)
ax.set_xlim(0, 27)
ax.set_ylim(-0.7, len(rows) + 1.7)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("k50U − k50F, exact-chain F1 points (each arm at its selected step 3072)")
style.legend_row(fig, 1000, 8, [("s", style.BLUE_LIGHT, style.BLUE_LIGHT, "seed 1"), ("s", style.BLUE, style.BLUE, "seed 2"),
                                 ("D", style.INK, style.INK, "two-seed mean")])
style.save(fig, "fig_p4seeds")
print("ok", pooled)
