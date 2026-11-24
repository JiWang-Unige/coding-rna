"""fig_ckpt  (SLOT_WIDE 1664x560)  slide c2-ckpt

Claim C2 (caveat): the P4 margin depends on which checkpoint each arm uses.

Source: docs/review_20261016/20_M34_stage1_results.md §13, stop-validation table lines 219-222
(k50F select, k50F new-species, k50U select / new-species, T10 select / new-species), the
selection line (k50F picks 3072) and the sensitivity sentence on line 243 (+39.3 / +6.12 /
-3.95 / +3.24; "+0.33 (4/8)" if k50F had been taken at step 2304). One seed per arm.
"select" = Arabidopsis + rice internal-validation regions (pooled); "new species" = the 8 new
species' internal-validation regions (pooled). Exact-chain F1 (%).
"""
import sys

sys.path.insert(0, "/Users/jiwang/Desktop/coding-rna-reorg/docs/review_deck_20261124/scripts")
import style  # noqa: E402

R20 = "docs/review_20261016/20_M34_stage1_results.md"

# ------------------------------------------------------------------ data
steps, f_sel, f_new, u_sel, u_new, t_sel, t_new = [], [], [], [], [], [], []
for ln in range(219, 223):
    c = style.cells(R20, ln)
    steps.append(int(c[0]))
    f_sel.append(style._num(c[1].split("（")[0]))
    f_new.append(style._num(c[2].split("（")[0]))
    a, b = c[3].split("/")
    u_sel.append(style._num(a)); u_new.append(style._num(b))
    a, b = c[4].split("/")
    t_sel.append(style._num(a)); t_new.append(style._num(b))
assert steps == [768, 1536, 2304, 3072]
gaps = [u - f for u, f in zip(u_new, f_new)]
style.has(R20, 243, "+39.3 / +6.12 / **−3.95** / +3.24", "+0.33（4/8 同向）")
for g, ref in zip(gaps, [39.3, 6.12, -3.95, 3.24]):
    style.close(g, ref, 0.04)
style.close(u_new[3] - f_new[2], 0.33, 0.006)       # k50U@3072 - k50F@2304
style.has(R20, 215, "**3072 步**", "65.76 / 69.30 / **69.76**")

# ------------------------------------------------------------------ figure
W, H = style.SLOT_WIDE
fig, ax0 = style.figure(W, H)
ax0.remove()
SER = [("k50U (marginalise)", style.BLUE, "o", "-"),
       ("k50F (fill background)", style.ORANGE, "o", "-"),
       ("T10 (all genes labelled)", style.GRAY, "D", "--")]


def panel(ax, series, ylim, yticks, title):
    for (label, col, mk, ls), ys in zip(SER, series):
        ax.plot(steps, ys, color=col, ls=ls, lw=2.6, zorder=2)
        ax.scatter(steps, ys, marker=mk, s=90 if mk == "o" else 60, zorder=3,
                   facecolor=style.BG if col == style.ORANGE else col, edgecolor=col, linewidth=2.4)
    ax.set_xticks(steps)
    ax.set_xlim(560, 3280)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_xlabel("Training step (stop validations)", labelpad=6)
    ax.set_title(title, pad=12)
    ax.axvspan(2930, 3214, color=style.GRID, alpha=0.6, zorder=0, lw=0)


axl = style.axes_px(fig, 92, 92, 560, 380)
panel(axl, [u_sel, f_sel, t_sel], (50, 76), range(50, 76, 5),
      "Selection score: Arabidopsis + rice")
axl.set_ylabel("Exact-chain F1 (%)", labelpad=6)
axl.text(3072, 51.2, "selected\n(all three arms)", ha="center", va="bottom", fontsize=15,
         color=style.MUTED, linespacing=1.1)

axr = style.axes_px(fig, 800, 92, 840, 380)
panel(axr, [u_new, f_new, t_new], (14, 86), range(20, 76, 10),
      "New-species internal validation: 8 species pooled")
axr.spines["left"].set_bounds(14, 75)
axr.text(470, 81, "k50U − k50F", ha="right", va="center", fontsize=16, color=style.MUTED,
         clip_on=False)
for s, g in zip(steps, gaps):
    ax_text = (f"{g:+.2f}" if abs(g) < 10 else f"{g:+.1f}").replace("-", "−")
    axr.text(s, 81, ax_text, ha="center", va="center", fontsize=17,
             fontweight="bold", color=style.BLUE if g > 0 else style.ORANGE)
axr.annotate("k50F at step 2304 vs k50U at 3072:\n+0.33 (4 of 8 species) — P4 would fail",
             xy=(2304, f_new[2]), xytext=(1180, 26), fontsize=16, color=style.BODY,
             ha="left", va="center", linespacing=1.15,
             arrowprops=dict(arrowstyle="-", color=style.MUTED, lw=1.4, shrinkB=8))

style.legend_row(fig, 86, 4, [
    ("o", style.BLUE, style.BLUE, SER[0][0]),
    ("o", style.BG, style.ORANGE, SER[1][0]),
    ("D", style.GRAY, style.GRAY, SER[2][0]),
    ("s", style.GRID, style.GRID, "selected checkpoint"),
], gap=34, fontsize=16)
print(style.save(fig, "fig_ckpt"))
